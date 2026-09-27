module ICE_IO_Probes
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ICE_Parameters_m
  use ICE_Base_Types_m,     only: ICE_vector_3D_type
  use ICE_Config_Types_m,   only: obj_io_probes, obj_sim_param
  use ICE_Advanced_Types_m, only: ICE_domain_type

  implicit none
  private
  public :: Setup_Probes, Write_Probes_Data

  type :: real_ptr
    real(R8), pointer :: p
  end type real_ptr

  type :: obj_probe
    type(ICE_vector_3D_type)          :: location
    integer, dimension(4)             :: ilocation
    integer                           :: nvar
    character(len=clen), allocatable  :: names(:)
    type(real_ptr), allocatable       :: variables(:)
    real(R8)                          :: dtime
    integer                           :: ntime = 0
    integer                           :: diter
    integer                           :: unit
  contains
    private
    procedure, pass(self) :: Place
  end type obj_probe

  integer, public :: nprobes = 0
  type(obj_probe), allocatable, target :: probe(:)

contains


  subroutine Setup_Probes(domain)
    use ICE_Config_Types_m, only: obj_io_probes, obj_sim_param
    use IR_Precision,       only: str
    use strings,            only: parse
#ifdef USE_MPI
    use ICE_Mod_MPI, only: is_local_block
#endif
    implicit none
    type(ICE_domain_type), intent(in) :: domain
    integer             :: i, ii, error
    character(len=clen) :: string(16)

    if (.not. allocated(obj_io_probes)) return
    nprobes = size(obj_io_probes)
    if (nprobes == 0) return
    allocate(probe(1:nprobes))

    do i = 1, nprobes

      obj_io_probes(i)%file = 'OUTPUT/'//trim(obj_io_probes(i)%file)//'.txt'

      string = ''
      probe(i)%nvar = 0
      call parse(obj_io_probes(i)%varnames, ' ', string)
      do ii = 1, size(string)
        if (trim(string(ii)) /= '') probe(i)%nvar = probe(i)%nvar + 1
      end do
      allocate(probe(i)%names(1:probe(i)%nvar))
      do ii = 1, probe(i)%nvar
        probe(i)%names(ii) = trim(string(ii))
      end do

      probe(i)%dtime = obj_io_probes(i)%dtime
      probe(i)%diter = obj_io_probes(i)%diter

      probe(i)%location%c = obj_io_probes(i)%loc
      probe(i)%ilocation  = obj_io_probes(i)%iloc

      if (sum(obj_io_probes(i)%iloc) == 0) call probe(i)%Place(domain)

#ifdef USE_MPI
      if (.not. is_local_block(probe(i)%ilocation(1))) cycle
#endif

      call Assign_Variables(probe(i), domain)

      if (.not. obj_sim_param%newrun) then
        open(newunit=probe(i)%unit, file=trim(obj_io_probes(i)%file), &
             status='OLD', iostat=error)
        if (error /= 0) then
          obj_io_probes(i)%error_message = &
            '[ERROR] Restart requested but probe file not found: '//trim(obj_io_probes(i)%file)
          write(*,'(A)') trim(obj_io_probes(i)%error_message)
          return
        end if
        error = 0
        do while (error == 0); read(probe(i)%unit, *, iostat=error); end do
        backspace(probe(i)%unit)
      else
        open(newunit=probe(i)%unit, file=trim(obj_io_probes(i)%file), &
             status='REPLACE', iostat=error)
      end if

    end do

  end subroutine Setup_Probes


  subroutine Place(self, domain)
    implicit none
    class(obj_probe),                  intent(inout) :: self
    type(ICE_domain_type), intent(in)                :: domain
    integer :: i, j, k, b
    real(R8) :: d0, d

    d0 = huge(1._R8)
    do b = 1, domain%nb
      do k = 0, domain%blk(b)%dim(3)
      do j = 0, domain%blk(b)%dim(2)
      do i = 1, domain%blk(b)%dim(1)
        d = norm2(self%location%c - domain%blk(b)%node(i,j,k)%c)
        if (d < d0) then
          d0 = d
          self%ilocation = [b, i, j, k]
        end if
      end do; end do; end do
    end do

  end subroutine Place


  subroutine Assign_Variables(probe, domain)
    use IR_Precision, only: str
    use ICE_Global_m, only: ngroups, ncond, nbase
    implicit none
    type(obj_probe),       intent(inout), target :: probe
    type(ICE_domain_type), intent(in),    target :: domain
    integer :: v, p, b, i, j, k

    b = probe%ilocation(1)
    i = probe%ilocation(2)
    j = probe%ilocation(3)
    k = probe%ilocation(4)

    allocate(probe%variables(1:probe%nvar))

    do v = 1, probe%nvar
      do p = 1, ngroups
        if (probe%names(v) == 'rho_'//trim(str(.true.,p))) &
          probe%variables(v)%p => domain%blk(b)%cond_phase(p)%prim(1,i,j,k)
        if (probe%names(v) == 'u_'//trim(str(.true.,p))) &
          probe%variables(v)%p => domain%blk(b)%cond_phase(p)%prim(2,i,j,k)
        if (probe%names(v) == 'v_'//trim(str(.true.,p))) &
          probe%variables(v)%p => domain%blk(b)%cond_phase(p)%prim(3,i,j,k)
        if (probe%names(v) == 'w_'//trim(str(.true.,p))) &
          probe%variables(v)%p => domain%blk(b)%cond_phase(p)%prim(4,i,j,k)
        if (probe%names(v) == 'T_'//trim(str(.true.,p))) &
          probe%variables(v)%p => domain%blk(b)%cond_phase(p)%prim(nbase(p)-1,i,j,k)
        if (probe%names(v) == 'n_'//trim(str(.true.,p))) &
          probe%variables(v)%p => domain%blk(b)%cond_phase(p)%prim(nbase(p),i,j,k)
      end do
    end do

  end subroutine Assign_Variables


  subroutine Write_Probes_Data(iter, time)
#ifdef USE_MPI
    use ICE_Mod_MPI, only: is_local_block
#endif
    implicit none
    integer, intent(in) :: iter
    real(R8), intent(in) :: time
    integer :: i, v

    do i = 1, nprobes
#ifdef USE_MPI
      if (.not. is_local_block(probe(i)%ilocation(1))) cycle
#endif
      if (mod(iter, probe(i)%diter) == 0) then
        write(probe(i)%unit,*) iter, (probe(i)%variables(v)%p, v=1,probe(i)%nvar)
      else if (time >= probe(i)%dtime * probe(i)%ntime) then
        write(probe(i)%unit,*) time, (probe(i)%variables(v)%p, v=1,probe(i)%nvar)
        probe(i)%ntime = probe(i)%ntime + 1
      end if
    end do

  end subroutine Write_Probes_Data


end module ICE_IO_Probes
