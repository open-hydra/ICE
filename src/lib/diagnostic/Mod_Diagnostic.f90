module ICE_Mod_Diagnostic
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ICE_Parameters_m

  implicit none
  private
  public :: Compute_Diagnostic, Write_Diagnostic

  character(len=llen), private :: Dvarnames = '"rho_p" "rhou_p" "rhov_p" "rhow_p" "rhoe_p" "dt"'

contains


  subroutine Compute_Diagnostic(new, old, dt, n, nc, average, total)
    use ICE_Global_m, only: gc, nres
    implicit none
    integer,  intent(in)    :: n(3), nc
    real(R8), intent(in)    :: new(nc, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc)
    real(R8), intent(in)    :: old(nc, 1-gc:n(1)+gc, 1-gc:n(2)+gc, 1-gc:n(3)+gc)
    real(R8), intent(in)    :: dt(n(1), n(2), n(3))
    real(R8), intent(out)   :: average(nres)
    real(R8), intent(inout) :: total(nres)
    ! Local
    integer  :: i, j, k
    real(R8) :: resn(nres), residuo(nres)

    residuo = 0._R8

    !$omp parallel
    !$omp do private(resn) reduction(+:residuo) collapse(3)
    do k = 1, n(3)
    do j = 1, n(2)
    do i = 1, n(1)
      resn(1)   = abs(new(1,    i,j,k) - old(1,    i,j,k))
      resn(2:4) = abs(new(2:4,  i,j,k) - old(2:4,  i,j,k))
      resn(5)   = abs(new(nc-1, i,j,k) - old(nc-1, i,j,k))
      residuo   = residuo + resn*resn
    end do; end do; end do
    !$omp end parallel

    average = sqrt(residuo / real(n(1)*n(2)*n(3), R8))
    total   = total + residuo

  end subroutine Compute_Diagnostic


  subroutine Write_Diagnostic(domain, IOfield, file)
    use ICE_Advanced_Types_m
    use ICE_Global_m
    use ICE_Config_Types_m,     only: obj_io
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
                 varnames=Dvarnames,                                     &
                 filename=trim(path)//trim(file)//'.tec')
      end select

    end if

    call mpi_io_barrier()

  end subroutine Write_Diagnostic


end module ICE_Mod_Diagnostic
