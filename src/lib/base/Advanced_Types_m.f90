module ICE_Advanced_Types_m
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ICE_Base_Types_m
  use ICE_Parameters_m
  use ICE_Series_Data_m
  use lib_orion_data, only: ORION_data

  implicit none

  !! ------------------------------------------------------
  !! FUNDAMENTAL TYPES
  !! ------------------------------------------------------
  type :: block_type
    integer                                        :: dim(3)   ! Cells in i-j-k (ghost excluded)
    real(R8),                          allocatable :: vol(:,:,:)
    type(ICE_vector_3D_type),          allocatable :: node(:,:,:)
    type(ICE_tensor_3D_type),          allocatable :: M(:,:,:)
    type(ICE_vector_3D_type),          allocatable :: dl(:,:,:)
    type(ICE_d_metrics_type)                       :: dir(3)
  end type block_type

  type :: bc_type
    integer  :: i, j, k, b, f                          ! Location in grid
    integer  :: type                                    ! BC type: ATLAS code (see ICE_IO_BC)
    integer  :: bs, is, js, ks, fs, d11, d12, d21, d22 ! Connection specs
    integer  :: partner = 0                             ! Connection: record of the same face seen from the other block
    logical  :: owner = .false.                         ! Connection: of the face's two records, the owner (Pair_Connections)
    integer  :: cf = 0                                  ! Connection: the face's number among its family's (Pair_Connections)
    integer               :: ni(2) = 0                 ! Chimera: donors per ghost layer
    integer,  allocatable :: donorID(:,:)               ! Chimera: (donor, [b i j k]), layer 1 then 2
    real(R8), allocatable :: volume_fraction(:)         ! Chimera: donor weights, sum to 1 per layer
    type(ICE_tensor_3D_type) :: Mg(2)                  ! Ghost cell metric tensor
    type(ICE_vector_3D_type) :: dlg(2)                 ! Ghost cell average cell length
    real(R8) :: volg(2)                                 ! Ghost cell volume
  end type bc_type

  !! ------------------------------------------------------
  !! ICE-SPECIFIC EXTENSIONS
  !! ------------------------------------------------------
  type :: ICE_gas_phase_type
    real(R8), dimension(:,:,:),   allocatable :: dt
    real(R8), dimension(:,:,:),   allocatable :: R, gam
    real(R8), dimension(:,:,:),   allocatable :: k, mu
    real(R8), dimension(:,:,:,:), allocatable :: prim
  end type ICE_gas_phase_type

  type :: ICE_cond_phase_type
    real(R8), dimension(:,:,:),     allocatable :: dt, tau, beta
    real(R8), dimension(:,:,:,:),   allocatable :: prim, prim_old
    real(R8), dimension(:,:,:,:),   allocatable :: source, residual
    ! First ghost cells that hold an injected state (ghost_record, codes 401-403), which the
    ! reconstructions place at the boundary face; over prim's cells, never set on an interior one
    logical,  dimension(:,:,:),     allocatable :: inj
  end type ICE_cond_phase_type

  type, extends(block_type) :: ICE_block_type
    type(ICE_gas_phase_type),                 allocatable :: gas_phase
    type(ICE_cond_phase_type), dimension(:),  allocatable :: cond_phase
  end type ICE_block_type

  !> The boundary cells of one family that this rank owns, each with its records
  !> in face order: cell c holds rec(first(c):first(c+1)-1).
  !> Built by build_local_bc_index; read by compute_bound.
  type :: ICE_bc_cells_type
    integer                            :: n = 0
    integer, dimension(:), allocatable :: first, rec
    ! With thread groups: the cells by owning group, tg_cell(tg_first(g):tg_first(g+1)-1) for group g
    integer, dimension(:), allocatable :: tg_cell, tg_first
  end type ICE_bc_cells_type

  !> The connection faces (101/201) of one family on one grid level. The flux of face cf is computed
  !> once, by its owner record, into F(:,cf): what enters the owner's cell through the face; the
  !> partner's cell takes its negative. Every rank numbers the faces alike (Pair_Connections), so a
  !> face's flux crosses ranks by its number.
  type :: ICE_conn_faces_type
    integer                               :: n = 0      ! faces of the family
    real(R8), dimension(:,:), allocatable :: F          ! (ncond, n)
    integer,  dimension(:),   allocatable :: prt        ! (n): the partner record of face cf
    real(R8), dimension(:),   allocatable :: bp         ! (n): the partner cell's detector weight, where that
                                                        !  cell lies on another rank (exchange_conn_beta)
    integer                               :: nown = 0   ! owner records in this rank's blocks
    integer,  dimension(:),   allocatable :: own        ! those records, in table order
    ! With thread groups: own by the group owning the record's block, own(tg_first(g):tg_first(g+1)-1) for group g
    integer,  dimension(:),   allocatable :: tg_first
  end type ICE_conn_faces_type

  type, extends(bc_type) :: ICE_bc_type
    integer  :: p                                       ! Particle group index
    real(R8) :: massflux, velocity, temperature
    real(R8) :: alpha, beta, radius
    real(R8) :: krho_sum = 0._R8                        ! 401: the krho of every family at the cell, summed
    type(time_series_type) :: BCtime
  end type ICE_bc_type

  !! ------------------------------------------------------
  !! COMPOUND TYPES
  !! ------------------------------------------------------
  type :: ICE_domain_type
    real(R8) :: time     = 0._R8
    real(R8) :: dtglobal = 0._R8
    integer  :: iter     = 0
    integer  :: itermax  = 0
    integer  :: nb       = 0
    integer  :: nbound   = 0
    ! Grid arrays
    integer, dimension(:,:), allocatable             :: n_bf
    type(ICE_block_type), dimension(:), allocatable  :: blk
    type(ICE_bc_type),    dimension(:), allocatable  :: bc
    ! MPI local BC indices (built by build_local_bc_index)
    integer                              :: n_local_bc = 0
    integer, dimension(:), allocatable  :: local_bc_idx
    integer                              :: n_local_bs = 0
    integer, dimension(:), allocatable  :: local_bs_idx
    type(ICE_bc_cells_type), dimension(:), allocatable :: bcells   ! one per family
    ! The same records grouped by family: family p is local_bc_grp(grp_first(p):grp_first(p+1)-1)
    integer, dimension(:), allocatable  :: local_bc_grp, grp_first
    ! With thread groups: positions in local_bc_grp by family and owning group -- family p, group g is
    ! local_bc_grp(tg_bc(tg_bc_first(g,p):tg_bc_first(g+1,p)-1))
    integer, dimension(:),   allocatable :: tg_bc
    integer, dimension(:,:), allocatable :: tg_bc_first
    type(ICE_conn_faces_type), dimension(:), allocatable :: conn   ! one per family: the connection faces' fluxes
  end type ICE_domain_type

  type :: ICE_simulation_type
    type(ICE_domain_type), allocatable :: domain(:)
    type(ORION_data)                   :: OCP
    type(ORION_data),      allocatable :: ODP(:)
  end type ICE_simulation_type

end module ICE_Advanced_Types_m
