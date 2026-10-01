module ICE_Global_m
  use ICE_Parameters_m
  implicit none

  character(len=clen) :: ICE_phase_prefix = 'part-'

  integer :: ndir     ! Number of spatial dimensions (set by setup_metrics)
  integer :: gc=2     ! Ghost cells per face
  integer :: nres=5   ! Number of tracked residuals (rho, u, v, w, T)
  real(R8), parameter :: rho_empty = 1e-6_R8   ! Bulk density below which a cell counts as empty
  integer :: tiles_per_thread = 1   ! Flux tiles per thread (ICE_TILES_PER_THREAD; a dynamic schedule needs more than one)

  ! Condensed-phase topology (set by Assign_Setup after reading input.ini)
  integer                                          :: ngroups
  integer                                          :: nrk
  integer,              dimension(:), allocatable  :: npop
  integer,              dimension(:), allocatable  :: ncond
  integer,              dimension(:), allocatable  :: nbase      ! slots of the closure alone
  logical,              dimension(:), allocatable  :: solid_of   ! family p's material solidifies
  integer, parameter                               :: ncond_max = 14   ! widest family: AG plus two slots

  ! Condensed materials (set by Setup_Materials): family p uses material mat_of(p)
  integer                                          :: nmat = 1
  integer,              dimension(:), allocatable  :: mat_of

end module ICE_Global_m
