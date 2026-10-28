!> Particles that enter the domain through extrapolation faces (400). The ghost of a 400 face copies its
!> boundary cell, so the face takes the flux of the cell's state whichever way the particles move: where they
!> move into the domain they enter with the state of that cell, and nothing outside the domain sets it. That
!> is right for the inlet of a truncated domain and wrong where a cloud is drawn in through an outlet, and the
!> state of the cell cannot tell the two apart: the run says where it happens. The flux is not changed.
!>
!> A boundary is the 400 faces of one face of a block, for one family, on one grid level. After the ghost fill
!> of each stage (one Runge-Kutta stage of one family's step), count_400_inflow counts the stages in which a
!> boundary cell that carries the phase (density above rho_empty) moves into the domain, w = v.n < 0 with n the
!> face normal out of the domain (the test of the inlets, ghost_rule), faster than round-off, and measures it:
!> |w| for MK, |w|/sqrt(P/rho) for IG and |w|/sqrt(Pnn/rho) for AG, Pnn = n.P.n. Faster than round-off is a
!> measure above 1e-12 for IG and AG, against their signal speed, so that a cell at rest whose whole velocity
!> is drift does not count, and -w above 1e-12 of the cell's speed for MK, which has no signal speed; for every
!> closure -w is also above the velocity of the closures' density floor. The rank that
!> owns the block prints and flushes a line at the first such stage of a boundary, so that a run that then stops
!> has it, and another each time the largest measure of a stage is more than 10 times that of the boundary's last
!> line, so that an entry of round-off does not hide a later inflow; report_400_inflow prints every boundary's
!> totals at the end of a standalone run, from root. A coupled run (hydra) never runs ICE's end-of-run block, so
!> there only the lines of the stages print. A restarted run counts from the restart.
module ICE_Mod_Inflow400
  use iso_fortran_env, only: R8 => real64, output_unit
  implicit none
  private
  public :: count_400_inflow, report_400_inflow

  !> The velocity the closures give a cell at their density floor, in each component (cons_2_prim_*): slower
  !> than that, only round-off moves the particles of a cell, such as the floor leaves next to an empty cell
  real(R8), parameter :: eps = 1e-25
  !> The part of the signal speed sqrt(P/rho) (IG, AG), or of the cell's speed (MK, which has no signal speed),
  !> below which a velocity into the domain is round-off
  real(R8), parameter :: round_off = 1e-12_R8

  type :: records_type
    integer, allocatable :: r(:)
  end type records_type

  !> The local 400 records of family p on level m, listed(p,m)%r: built at the first stage of level m
  type(records_type), allocatable :: listed(:,:)
  !> Per boundary (face, block, family, level): the stages with particles entering, the first one's iteration
  !> and faces, the largest measure at that stage and over the run, and the largest measure of the boundary's
  !> last line (zero before its first)
  integer,  allocatable :: nstages(:,:,:,:), first_iter(:,:,:,:), first_faces(:,:,:,:)
  real(R8), allocatable :: first_x(:,:,:,:), largest_x(:,:,:,:), line_x(:,:,:,:)

contains


  !> The stage of family p on the current level, after its ghost fill: its local 400 records, in one serial
  !> pass (a rank's records are those of its own blocks, so its counts are the whole of its boundaries'),
  !> folded into the boundaries' counts and lines. Called from serial code, between the parallel regions of the
  !> stage: a line prints once per boundary and stage, from the counts over all its faces, whatever the threads.
  subroutine count_400_inflow(grid, p)
    use ICE_Advanced_Types_m, only: ICE_domain_type
    use ICE_Global_m,         only: rho_empty
    use ICE_Config_Types_m,   only: obj_multigrid, obj_sim_param, obj_time_scheme
    implicit none
    type(ICE_domain_type), intent(in) :: grid
    integer,               intent(in) :: p
    integer  :: lev, ii, n, b, f, i, j, k, nface(6, grid%nb), nabove(6, grid%nb)
    real(R8) :: q(10), normal(3), w, x, pnn, xs(6, grid%nb), xa(6, grid%nb)
    character(len=48) :: which

    lev = obj_multigrid%MG_level
    call list_records(grid, lev)
    if (size(listed(p,lev)%r) == 0) return

    nface = 0 ; nabove = 0
    xs    = 0._R8 ; xa = 0._R8
    do ii = 1, size(listed(p,lev)%r)
      n = listed(p,lev)%r(ii)
      b = grid%bc(n)%b ; i = grid%bc(n)%i ; j = grid%bc(n)%j ; k = grid%bc(n)%k ; f = grid%bc(n)%f
      q(1:4) = grid%blk(b)%cond_phase(p)%prim(1:4,i,j,k)
      if (.not. (q(1) > rho_empty)) cycle
      if (f <= 2) then
        normal = grid%blk(b)%dir(1)%f(i-mod(f,2),j,k)%N
      elseif (f <= 4) then
        normal = grid%blk(b)%dir(2)%f(i,j-mod(f,2),k)%N
      else
        normal = grid%blk(b)%dir(3)%f(i,j,k-mod(f,2))%N
      endif
      normal = real(1-2*mod(f,2), R8) * normal
      w = dot_product(q(2:4), normal)
      if (.not. (w < -eps)) cycle
      select case (trim(obj_time_scheme%model(p)))
      case ('IG')
        q(5) = grid%blk(b)%cond_phase(p)%prim(5,i,j,k)
        x = -w * sqrt(q(1) / q(5))
        if (.not. (x > round_off)) cycle
      case ('AG')
        q(5:10) = grid%blk(b)%cond_phase(p)%prim(5:10,i,j,k)
        pnn = q(5)*normal(1)**2 + q(8)*normal(2)**2 + q(10)*normal(3)**2 &
            + 2._R8*(q(6)*normal(1)*normal(2) + q(7)*normal(1)*normal(3) + q(9)*normal(2)*normal(3))
        x = -w * sqrt(q(1) / max(pnn, tiny(1._R8)))
        if (.not. (x > round_off)) cycle
      case default
        x = -w
        if (.not. (x > round_off*norm2(q(2:4)))) cycle
      end select
      nface(f,b) = nface(f,b) + 1
      xs(f,b)    = max(xs(f,b), x)
      if (x > 10._R8*line_x(f,b,p,lev)) then
        nabove(f,b) = nabove(f,b) + 1
        xa(f,b)     = max(xa(f,b), x)
      endif
    enddo

    do b = 1, grid%nb
      do f = 1, 6
        if (nface(f,b) == 0) cycle
        if (nstages(f,b,p,lev) == 0) then
          first_iter(f,b,p,lev)  = obj_sim_param%iter_general
          first_faces(f,b,p,lev) = nface(f,b)
          first_x(f,b,p,lev)     = xs(f,b)
        endif
        nstages(f,b,p,lev)   = nstages(f,b,p,lev) + 1
        largest_x(f,b,p,lev) = max(largest_x(f,b,p,lev), xs(f,b))
        if (nabove(f,b) == 0) cycle
        which = ', faces: '
        if (line_x(f,b,p,lev) > 0._R8) which = ', faces above 10 times its last line: '
        write(*,'(A,I0,A,I0,A)') ' 400 boundary of '//trim(boundary_name(b, f, p, lev))// &
          ': particles entering at iteration ', obj_sim_param%iter_general, trim(which)//' ', nabove(f,b), &
          ', largest '//trim(measure_text(p, xa(f,b)))
        flush(output_unit)
        line_x(f,b,p,lev) = xa(f,b)
      enddo
    enddo

  end subroutine count_400_inflow


  !> The boundaries through which particles entered over the run, one line each: the stages, the first one's
  !> iteration, faces and largest measure, and the largest measure over the run; nothing when none entered.
  !> Every rank calls it, in every run: one reduction to root, which prints. A boundary is counted only by the
  !> rank that owns its block and is zero on the others, so the sum over the ranks is that rank's numbers, and
  !> a parallel run prints the totals of a serial one.
  subroutine report_400_inflow(nb)
    use ICE_Global_m,       only: ngroups
    use ICE_Config_Types_m, only: obj_multigrid
    use ICE_Mod_MPI,        only: mpi_is_root, mpi_reduce_sum_r8_array
    implicit none
    integer, intent(in) :: nb
    real(R8), allocatable :: buf(:,:)
    integer :: nk, b, f, p, lev, m

    if (.not. allocated(nstages)) call allocate_counts(nb, ngroups, obj_multigrid%MGL)
    nk = size(nstages)
    allocate(buf(nk, 5))
    buf(:,1) = real(reshape(nstages,     [nk]), R8)
    buf(:,2) = real(reshape(first_iter,  [nk]), R8)
    buf(:,3) = real(reshape(first_faces, [nk]), R8)
    buf(:,4) = reshape(first_x,   [nk])
    buf(:,5) = reshape(largest_x, [nk])
    call mpi_reduce_sum_r8_array(buf, 5*nk)
    if (.not. mpi_is_root) return

    do p = 1, size(nstages, 3)
      do lev = size(nstages, 4), 1, -1
        do b = 1, size(nstages, 2)
          do f = 1, 6
            m = f + 6*((b-1) + size(nstages, 2)*((p-1) + size(nstages, 3)*(lev-1)))
            if (buf(m,1) == 0._R8) cycle
            write(*,'(A,I0,A,I0,A,I0,A)') ' 400 boundary of '//trim(boundary_name(b, f, p, lev))// &
              ': particles entered in ', nint(buf(m,1)), ' stages, first at iteration ', nint(buf(m,2)), &
              ' (faces: ', nint(buf(m,3)), ', largest '//trim(measure_text(p, buf(m,4)))// &
              '), largest over the run: '//trim(number_text(p, buf(m,5)))
          enddo
        enddo
      enddo
    enddo

  end subroutine report_400_inflow


  !> The counts, zero, for nb blocks, ng families and nl levels
  subroutine allocate_counts(nb, ng, nl)
    implicit none
    integer, intent(in) :: nb, ng, nl
    allocate(nstages(6,nb,ng,nl), first_iter(6,nb,ng,nl), first_faces(6,nb,ng,nl), source=0)
    allocate(first_x(6,nb,ng,nl), largest_x(6,nb,ng,nl), line_x(6,nb,ng,nl), source=0._R8)
    allocate(listed(ng,nl))
  end subroutine allocate_counts


  !> At the first stage of level lev, the local 400 records of each family on it
  subroutine list_records(grid, lev)
    use ICE_Advanced_Types_m, only: ICE_domain_type
    use ICE_Global_m,         only: ngroups
    use ICE_Config_Types_m,   only: obj_multigrid
    implicit none
    type(ICE_domain_type), intent(in) :: grid
    integer,               intent(in) :: lev
    integer :: p, ii, lo, hi

    if (.not. allocated(nstages)) call allocate_counts(grid%nb, ngroups, obj_multigrid%MGL)
    if (allocated(listed(1,lev)%r)) return
    do p = 1, ngroups
      lo = grid%grp_first(p) ; hi = grid%grp_first(p+1) - 1
      listed(p,lev)%r = pack(grid%local_bc_grp(lo:hi), [(grid%bc(grid%local_bc_grp(ii))%type == 400, ii = lo, hi)])
    enddo
  end subroutine list_records


  !> "block b, face f", with the family when the run has several and the level when it has several
  function boundary_name(b, f, p, lev) result(s)
    use ICE_Global_m,       only: ngroups
    use ICE_Config_Types_m, only: obj_multigrid
    implicit none
    integer, intent(in) :: b, f, p, lev
    character(len=80)   :: s
    write(s,'(A,I0,A,I0)') 'block ', b, ', face ', f
    if (ngroups > 1)           write(s,'(A,I0)') trim(s)//', family ', p
    if (obj_multigrid%MGL > 1) write(s,'(A,I0)') trim(s)//', grid level ', lev
  end function boundary_name


  !> "|w|: <x> m/s" (MK), "|w|/sqrt(P/rho): <x>" (IG) or "|w|/sqrt(Pnn/rho): <x>" (AG), for family p
  function measure_text(p, x) result(s)
    use ICE_Config_Types_m, only: obj_time_scheme
    implicit none
    integer,  intent(in) :: p
    real(R8), intent(in) :: x
    character(len=60)    :: s
    select case (trim(obj_time_scheme%model(p)))
    case ('IG')   ; s = '|w|/sqrt(P/rho): '//number_text(p, x)
    case ('AG')   ; s = '|w|/sqrt(Pnn/rho): '//number_text(p, x)
    case default  ; s = '|w|: '//number_text(p, x)
    end select
  end function measure_text


  !> x, with its unit for family p: m/s for MK, none for IG and AG. An exponent of three digits gets its E too
  !> (an AG pressure along the normal at zero gives a measure of about 1e150)
  function number_text(p, x) result(s)
    use ICE_Config_Types_m, only: obj_time_scheme
    implicit none
    integer,  intent(in) :: p
    real(R8), intent(in) :: x
    character(len=30)    :: s
    if (x >= 1e100_R8 .or. (x > 0._R8 .and. x < 1e-99_R8)) then
      write(s,'(ES12.4E3)') x
    else
      write(s,'(ES11.4)') x
    endif
    s = adjustl(s)
    select case (trim(obj_time_scheme%model(p)))
    case ('IG', 'AG')
    case default ; s = trim(s)//' m/s'
    end select
  end function number_text

end module ICE_Mod_Inflow400
