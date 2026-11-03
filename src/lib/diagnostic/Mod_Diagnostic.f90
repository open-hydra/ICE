module ICE_Mod_Diagnostic
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ICE_Parameters_m

  implicit none
  private
  public :: residual_norm, residual_columns, write_residual_row, Write_Diagnostic

contains


  !> The residual norm: the L2 norm of prim - prim_old in density, momentum and
  !> energy over every cell, block and family of the domain. The sum runs in one
  !> fixed order -- the cells of a block in k, j, i, then the blocks and families
  !> in their order -- on every rank: a rank sums its own blocks, the per-block
  !> partial sums are exchanged (one contributor per entry, so the exchange is
  !> exact), and every rank adds them up in the same order. The thread count
  !> and the rank count therefore do not change a bit of it, which a reduction
  !> over threads and a sum over ranks did. With family, the same norms of each
  !> family alone, from the same partial sums added over the blocks in the same
  !> order (total is not derived from them, so it keeps its bits).
  subroutine residual_norm(grid, total, family)
    use ICE_Global_m,           only: ngroups, ncond, nbase, nres
    use ICE_Advanced_Types_m,   only: ICE_domain_type
    use ICE_Mod_MPI,            only: is_local_block, mpi_allreduce_sum_r8_array
    use ICE_Mod_Timers,         only: timer_sync_begin, timer_sync_end
    implicit none
    type(ICE_domain_type), intent(in)            :: grid
    real(R8),              intent(out)           :: total(nres)
    real(R8),              intent(out), optional :: family(nres, ngroups)
    ! Local
    real(R8), allocatable :: part(:,:)
    integer :: b, p, m

    allocate(part(nres, grid%nb*ngroups))
    part = 0._R8
    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle
      do p = 1, ngroups
        m = (b-1)*ngroups + p
        call residual_sq(grid%blk(b)%cond_phase(p)%prim, grid%blk(b)%cond_phase(p)%prim_old, &
                         grid%blk(b)%dim, ncond(p), nbase(p)-1, part(:,m))
      end do
    end do

    call timer_sync_begin()
    call mpi_allreduce_sum_r8_array(part, nres*grid%nb*ngroups)
    call timer_sync_end()

    total = 0._R8
    do b = 1, grid%nb
      do p = 1, ngroups
        total = total + part(:, (b-1)*ngroups + p)
      end do
    end do
    total = sqrt(total)

    if (present(family)) then
      family = 0._R8
      do b = 1, grid%nb
        do p = 1, ngroups
          family(:, p) = family(:, p) + part(:, (b-1)*ngroups + p)
        end do
      end do
      family = sqrt(family)
    end if
    deallocate(part)

  end subroutine residual_norm


  !> The values of a residual-history.dat row after the iteration and the time:
  !> the norms of the domain, then, with several families, the same norms of
  !> each family in turn.
  integer function residual_columns()
    use ICE_Global_m, only: ngroups, nres
    implicit none

    residual_columns = nres
    if (ngroups > 1) residual_columns = nres * (1 + ngroups)

  end function residual_columns


  !> One row of residual-history.dat, flushed: the global iteration, the time and
  !> the norms (residual_columns). The root writes it.
  subroutine write_residual_row(unit, fmt, iter, time)
    use ICE_Global_m,       only: ngroups
    use ICE_Config_Types_m, only: obj_sim_param
    implicit none
    integer,          intent(in) :: unit, iter
    character(len=*), intent(in) :: fmt
    real(R8),         intent(in) :: time

    if (ngroups > 1) then
      write(unit, fmt) iter, time, obj_sim_param%residuotot, obj_sim_param%residuofam
    else
      write(unit, fmt) iter, time, obj_sim_param%residuotot
    end if
    flush(unit)

  end subroutine write_residual_row


  !> Sum of the squared changes of one family on one block, cells in k, j, i order.
  subroutine residual_sq(new, old, n, nc, iT, residuo)
    use ICE_Global_m, only: gc, nres
    implicit none
    integer,  intent(in)  :: n(3), nc, iT
    real(R8), intent(in)  :: new(nc, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc)
    real(R8), intent(in)  :: old(nc, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc)
    real(R8), intent(out) :: residuo(nres)
    ! Local
    integer  :: i, j, k
    real(R8) :: resn(nres)

    residuo = 0._R8
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)
      resn(1)   = abs(new(1,    i,j,k) - old(1,    i,j,k))
      resn(2:4) = abs(new(2:4,  i,j,k) - old(2:4,  i,j,k))
      resn(5)   = abs(new(iT, i,j,k) - old(iT, i,j,k))
      residuo   = residuo + resn*resn
    end do; end do; end do

  end subroutine residual_sq


  !> The residuals of rho, rho u, rho v, rho w and rho e and the local time step
  !> of every family, one named column each (rho_p<n> ... dt_p<n>), in the
  !> output buffer of the field.
  subroutine Write_Diagnostic(domain, IOfield, file)
    use ICE_Advanced_Types_m
    use ICE_Global_m
    use ICE_Config_Types_m,     only: obj_io
    use ICE_IO_Solution,        only: io_extension
    use ICE_Mod_MPI,            only: mpi_is_root, mpi_abort_all
    use ICE_Mod_GhostExchange,  only: gather_diagnostic_to_root, mpi_io_barrier
    use IR_Precision,           only: str
    use strings,                only: parse
    use Lib_VTK
    use Lib_Tecplot
    implicit none
    type(ICE_domain_type), intent(inout) :: domain
    type(orion_data),      intent(inout) :: IOfield
    character(len=*),      intent(in)    :: file
    ! Local
    character(len=llen) :: path, localpath_vtk, written
    character(len=clen) :: format(2)
    character(len=hlen) :: varnames, vtknames, what
    character(len=8)    :: f
    real(R8), allocatable :: full(:,:,:,:)
    integer(kind=I4)    :: E_IO, b, p, c, nstart, nend, ndiag, nfield

    call gather_diagnostic_to_root(domain)

    if (mpi_is_root) then

      path = 'OUTPUT/'
      call parse(obj_io%sol_format, ' ', format)

      ndiag = 6*ngroups
      varnames = ''
      do p = 1, ngroups
        f = '_p'//trim(str(.true.,p))
        varnames = trim(varnames)//' "rho'//trim(f)//'" "rhou'//trim(f)//'" "rhov'//trim(f)//'" "rhow'//trim(f)// &
                   '" "rhoe'//trim(f)//'" "dt'//trim(f)//'"'
      end do

      do b = 1, domain%nb
        nend = 0
        do p = 1, ngroups
          nstart = nend + 1;  nend = nstart - 1 + 6
          IOfield%block(b)%vars(nstart:nend-2,:,:,:) = &
            domain%blk(b)%cond_phase(p)%residual(1:4, 1:domain%blk(b)%dim(1), &
                                                     1:domain%blk(b)%dim(2), &
                                                     1:domain%blk(b)%dim(3))
          IOfield%block(b)%vars(nend-1, 1:IOfield%block(b)%Ni, &
                                        1:IOfield%block(b)%Nj, &
                                        1:IOfield%block(b)%Nk) = &
            domain%blk(b)%cond_phase(p)%residual(5, 1:domain%blk(b)%dim(1), &
                                                   1:domain%blk(b)%dim(2), &
                                                   1:domain%blk(b)%dim(3))
          IOfield%block(b)%vars(nend,   1:IOfield%block(b)%Ni, &
                                        1:IOfield%block(b)%Nj, &
                                        1:IOfield%block(b)%Nk) = &
            domain%blk(b)%cond_phase(p)%dt(1:domain%blk(b)%dim(1), &
                                          1:domain%blk(b)%dim(2), &
                                          1:domain%blk(b)%dim(3))
        end do
      end do

      E_IO = 0
      select case(trim(format(1)))
      case('vtk')
        IOfield%vtk%format = trim(format(2))
        ! The names become XML attributes (Name="..."): the quotes of the list go
        vtknames = varnames
        do c = 1, len_trim(vtknames)
          if (vtknames(c:c) == '"') vtknames(c:c) = ' '
        end do
        ! The vtk writer writes every column the buffer holds, and the buffer of
        ! the field has a column per field variable: it gets the diagnostic
        ! columns only, then the field's count back (the field writer fills
        ! every column at each write).
        nfield = size(IOfield%block(1)%vars, 1)
        call keep_columns(ndiag)
        localpath_vtk = trim(path)//'vtk/'
        call execute_command_line('mkdir -p '//trim(localpath_vtk))
        written = trim(path)//trim(file)//'.vtm'
        ! TIME as the solution's (write_vtk_tec)
        E_IO = vtk_write_structured_multiblock(orion=IOfield,            &
                 vtspath=trim(localpath_vtk)//trim(file),                &
                 vtmpath=trim(path)//trim(file),                         &
                 varnames=vtknames, time=IOfield%solutiontime)
        call keep_columns(nfield)
      case('tecplot')
        IOfield%tec%format = trim(format(2))
        written = trim(path)//trim(file)//io_extension(format)
        E_IO = tec_write_structured_multiblock(orion=IOfield,            &
                 varnames=varnames, Nvars=ndiag,                         &
                 filename=trim(written))
      end select

      ! A write that the library refuses or that fails stops the run, with the
      ! file (as the solution's, write_vtk_tec)
      if (E_IO /= 0) then
        write(what,'(A,I0,A)') 'cannot write '//trim(written)//' (error ', E_IO, ' of the writer)'
        call mpi_abort_all(trim(what))
      end if

    end if

    call mpi_io_barrier()

  contains

    !> Every block of the buffer with n columns, the first of them kept.
    subroutine keep_columns(n)
      integer, intent(in) :: n
      integer :: b, m

      do b = 1, size(IOfield%block)
        call move_alloc(IOfield%block(b)%vars, full)
        allocate(IOfield%block(b)%vars(n, lbound(full,2):ubound(full,2), lbound(full,3):ubound(full,3), &
                                          lbound(full,4):ubound(full,4)))
        m = min(n, size(full, 1))
        IOfield%block(b)%vars(1:m,:,:,:) = full(1:m,:,:,:)
        deallocate(full)
      end do
    end subroutine keep_columns

  end subroutine Write_Diagnostic


end module ICE_Mod_Diagnostic
