module ICE_Mod_Diagnostic
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ICE_Parameters_m

  implicit none
  private
  public :: residual_norm, Write_Diagnostic

  character(len=llen), private :: Dvarnames = '"rho_p" "rhou_p" "rhov_p" "rhow_p" "rhoe_p" "dt"'

contains


  !> The residual norm: the L2 norm of prim - prim_old in density, momentum and
  !> energy over every cell, block and family of the domain. The sum runs in one
  !> fixed order -- the cells of a block in k, j, i, then the blocks and families
  !> in their order -- on every rank: a rank sums its own blocks, the per-block
  !> partial sums are exchanged (one contributor per entry, so the exchange is
  !> exact), and every rank adds them up in the same order. The thread count
  !> and the rank count therefore do not change a bit of it, which a reduction
  !> over threads and a sum over ranks did.
  subroutine residual_norm(grid, total)
    use ICE_Global_m,           only: ngroups, ncond, nbase, nres
    use ICE_Advanced_Types_m,   only: ICE_domain_type
    use ICE_Mod_MPI,            only: is_local_block, mpi_allreduce_sum_r8_array
    use ICE_Mod_Timers,         only: timer_sync_begin, timer_sync_end
    implicit none
    type(ICE_domain_type), intent(in)  :: grid
    real(R8),              intent(out) :: total(nres)
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
    deallocate(part)

  end subroutine residual_norm


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


  subroutine Write_Diagnostic(domain, IOfield, file)
    use ICE_Advanced_Types_m
    use ICE_Global_m
    use ICE_Config_Types_m,     only: obj_io
    use ICE_IO_Solution,        only: io_extension
    use ICE_Mod_MPI,            only: mpi_is_root
    use ICE_Mod_GhostExchange,  only: gather_diagnostic_to_root, mpi_io_barrier
    use strings,                only: parse
    use Lib_VTK
    use Lib_Tecplot
    implicit none
    type(ICE_domain_type), intent(inout) :: domain
    type(orion_data),      intent(inout) :: IOfield
    character(len=*),      intent(in)    :: file
    ! Local
    character(len=llen) :: path, localpath_vtk
    character(len=clen) :: format(2)
    integer(kind=I4)    :: E_IO, b, p, nstart, nend

    call gather_diagnostic_to_root(domain)

    if (mpi_is_root) then

      path = 'OUTPUT/'
      call parse(obj_io%sol_format, ' ', format)

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

      select case(trim(format(1)))
      case('vtk')
        IOfield%vtk%format = trim(format(2))
        localpath_vtk = trim(path)//'vtk/'
        call execute_command_line('mkdir -p '//trim(localpath_vtk))
        E_IO = vtk_write_structured_multiblock(orion=IOfield,            &
                 vtspath=trim(localpath_vtk)//trim(file),                &
                 vtmpath=trim(path)//trim(file),                         &
                 varnames=Dvarnames, time=domain%time)
      case('tecplot')
        IOfield%tec%format = trim(format(2))
        E_IO = tec_write_structured_multiblock(orion=IOfield,            &
                 varnames=Dvarnames, Nvars=6*ngroups,                    &
                 filename=trim(path)//trim(file)//io_extension(format))
      end select

    end if

    call mpi_io_barrier()

  end subroutine Write_Diagnostic


end module ICE_Mod_Diagnostic
