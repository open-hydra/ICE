module ICE_Mod_Multigrid
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use ICE_Mod_MPI, only: is_local_block

  implicit none
  private
  public :: Setup_Multigrid, Restriction, Prolongation

contains


  !> Build every coarse level and leave it ready to be solved: its own grid, metrics,
  !> boundary table, MPI halo schedule and a restriction of the initial solution.
  !> Ownership is decided once, on the fine level, by partition_blocks -- a coarse level
  !> has the same blocks in the same order, so the same rank owns the same block on all
  !> of them and no second partition is needed. The halo schedule is NOT shared: it names
  !> cells, and a coarse cell is not a fine one, so each level builds its own.
  subroutine Setup_Multigrid(simulation)
    use ICE_Advanced_Types_m
    use ICE_Config_Types_m,      only: obj_multigrid
    use ICE_Global_m,            only: ngroups
    use ICE_IO_BC,               only: Setup_BC
    use ICE_Lib_Ghost,           only: compute_ghost, fill_second_ghost
    use ICE_Lib_Multigrid,       only: Check_Multigrid, Coarse_Grid, Coarse_IOfield
    use ICE_Mod_Allocate_Data,   only: Allocate_Block
    use ICE_Mod_GhostExchange,   only: build_ghost_schedule, build_local_bc_index, &
                                       select_ghost_level
    use ICE_Mod_Metrics,         only: setup_metrics
    implicit none
    type(ICE_simulation_type), intent(inout) :: simulation
    integer :: b, m, i, j, k, nb, rap, nbound, ni, nj, nk

    call Check_Multigrid(simulation%domain(1))

    nb = simulation%domain(1)%nb

    do m = 2, obj_multigrid%MGL

      allocate(simulation%domain(m)%blk(nb))
      simulation%domain(m)%nb   = nb
      simulation%domain(m)%time = simulation%domain(1)%time

      rap = 2
      nbound = 0
      do b = 1, nb
        simulation%domain(m)%blk(b)%dim = simulation%domain(m-1)%blk(b)%dim / rap
        simulation%domain(m)%blk(b)%dim(3) = max(1, simulation%domain(m)%blk(b)%dim(3))
        call Allocate_Block(simulation%domain(m)%blk(b), simulation%domain(m)%blk(b)%dim, b)

        ni = simulation%domain(m)%blk(b)%dim(1)
        nj = simulation%domain(m)%blk(b)%dim(2)
        nk = simulation%domain(m)%blk(b)%dim(3)
        nbound = nbound + 2*nj*nk + 2*ni*nk + 2*nj*ni
      end do

      allocate(simulation%domain(m)%bc(1:nbound*ngroups))
      allocate(simulation%domain(m)%n_bf(1:nb, 1:6))

      call Coarse_Grid(simulation%domain(m-1), simulation%domain(m))

      call Coarse_IOfield(simulation%ODP(m-1), simulation%ODP(m))

      do b = 1, nb
        do k = 0, simulation%domain(m)%blk(b)%dim(3)
        do j = 0, simulation%domain(m)%blk(b)%dim(2)
        do i = 0, simulation%domain(m)%blk(b)%dim(1)
          simulation%ODP(m)%block(b)%mesh(:,i,j,k) = simulation%domain(m)%blk(b)%node(i,j,k)%c
        end do ; end do ; end do
      end do

      call setup_metrics(simulation%domain(m), simulation%ODP(m))

      call Setup_BC(simulation%domain(m), m)
      call build_local_bc_index(simulation%domain(m))
      call build_ghost_schedule(simulation%domain(m), m)

      call Restriction(simulation%domain(m-1), simulation%domain(m))
      call select_ghost_level(m)
      call compute_ghost(simulation%domain(m))
      call fill_second_ghost(simulation%domain(m))

    end do

    ! Set iteration limits per level
    do m = 1, obj_multigrid%MGL
      simulation%domain(m)%itermax = obj_multigrid%iter_threshold(m)
    end do

  end subroutine Setup_Multigrid


  subroutine Restriction(Fine, Coarse)
    use ICE_Advanced_Types_m
    use ICE_Global_m,      only: ngroups
    use ICE_Lib_Multigrid, only: fine2coarse_prim
    implicit none
    type(ICE_domain_type), intent(inout) :: Fine, Coarse
    integer :: b, p

    do b = 1, Coarse%nb
      do p = 1, ngroups
        call fine2coarse_prim(p,                                     &
          Fine%blk(b)%cond_phase(p)%prim,                            &
          Coarse%blk(b)%cond_phase(p)%prim,                          &
          Fine%blk(b)%vol,   Coarse%blk(b)%vol,                      &
          Fine%blk(b)%dim,   Coarse%blk(b)%dim)
      end do
    end do

  end subroutine Restriction


  subroutine Prolongation(Fine, Coarse)
    use ICE_Advanced_Types_m
    use ICE_Global_m,      only: ngroups
    use ICE_Lib_Multigrid, only: coarse2fine_prim
    implicit none
    type(ICE_domain_type), intent(inout) :: Fine, Coarse
    integer :: b, p

    do b = 1, Coarse%nb
      if (.not. is_local_block(b)) cycle
      do p = 1, ngroups
        call coarse2fine_prim(p,                                     &
          Fine%blk(b)%cond_phase(p)%prim,                            &
          Coarse%blk(b)%cond_phase(p)%prim,                          &
          Fine%blk(b)%dim,   Coarse%blk(b)%dim)
      end do
    end do

  end subroutine Prolongation


end module ICE_Mod_Multigrid
