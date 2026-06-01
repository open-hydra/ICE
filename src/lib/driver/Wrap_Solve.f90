module ICE_Wrap_Solve
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: ICE_solve

contains

  subroutine ICE_solve(simulation, External_Function)
    use ICE_Advanced_Types_m, only: ICE_simulation_type
    use ICE_Config_Types_m,   only: obj_multigrid, obj_sim_param
    use ICE_Mod_Explicit,     only: Explicit_Step
    use ICE_Mod_Multigrid,    only: Prolongation
    implicit none
    type(ICE_simulation_type), intent(inout) :: simulation
    external :: External_Function
    integer :: lv

    lv = obj_multigrid%MG_level
    call Explicit_Step(simulation%domain(lv))

    ! Coarse level exhausted its iteration budget → trigger level transition
    if (lv > 1 .and. obj_sim_param%iter_from_call >= simulation%domain(lv)%itermax) then
      obj_multigrid%change_MG = .true.
      ! Don't decrement MG_level here — done in Wrap_Postprocess (matching MOSE)
      if (obj_sim_param%TODO == 3) obj_sim_param%TODO = 2
    end if

    ! Prolongation onto finer grid (change_MG reset in Wrap_Postprocess)
    if (obj_multigrid%change_MG) then
      call Prolongation(Fine   = simulation%domain(lv-1), &
                        Coarse = simulation%domain(lv))
    end if

  end subroutine ICE_solve

end module ICE_Wrap_Solve
