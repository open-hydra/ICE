!> ICE's exchange with the solvers coupled to it: the only module a coupler needs.
!>
!> What crosses ICE's boundary goes by one of four routes, cell by cell:
!>
!>   gas state   gas -> ICE     set_gas            the gas that the drops of a cell see
!>   exchange    ICE -> gas     exchange_rate      what a family gives the gas, per unit time, over its last step
!>   removals    ICE -> gas     take_removed       what the burn-out removed from a cell at once: an amount
!>   births      any -> ICE     add_birth          a mass born at once in a family, for instance another solver's
!>                              add_birth_rate     a rate of births, in force until clear_birth_rates
!>
!> Two calls at set-up, after ICE's own set-up and before its first step, declare where mass is born:
!>
!>   set_generation   a burning material generates a condensed mass, born as drops of another family
!>   birth_family     the family in which another solver's mass of a given material is born
!>
!> ICE's own generated mass takes the births route too (ICE_Mod_Sources), so a mass is born in the same way whoever
!> generates it.
!>
!> Units. Amounts are per unit volume (kg/m^3, kg/(m^2 s), J/m^3) and rates per unit volume and time. Energies are
!> ICE's: the energy e(T) of the material on ICE's datum, plus the kinetic energy, as the conservative state holds
!> them. A coupler whose gas measures energy from another datum adds its own offsets.
!>
!> A birth is described by a packet of six numbers: [mass, momentum (3), energy, number of drops] (birth_packet).
module ICE_Mod_Exchange
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private

  public :: exchange_t
  public :: set_gas, exchange_rate, take_removed
  public :: set_generation, birth_family
  public :: birth_packet, add_birth, add_birth_rate, clear_birth_rates

  !> What one family gives out of one cell: per unit volume and time for exchange_rate, per unit volume for
  !> take_removed. The drops lose mass by evaporating or burning. Part of it can be born again as drops of another
  !> family, the generated mass (the oxide of a burning metal, for instance); the gas receives the rest.
  type :: exchange_t
    real(R8) :: mass        = 0._R8   ! the mass the gas receives
    real(R8) :: momentum(3) = 0._R8   ! the momentum the gas receives: the drag, and the momentum of that mass
    real(R8) :: energy      = 0._R8   ! the energy the gas receives: the heat exchanged, and the energy of that mass
    real(R8) :: generated   = 0._R8   ! the mass born as drops of the generation family instead of reaching the gas
    real(R8) :: heat        = 0._R8   ! the heat of reaction of the mass that burns (see below)
  end type exchange_t
  ! heat: a gas made of one species adds it to its energy. A gas whose species carry their formation enthalpies
  ! already gains it from the change of species, and leaves it out.

contains

  ! ---------------------------------------------------------------------------------------------------------------
  ! Gas and ICE
  ! ---------------------------------------------------------------------------------------------------------------

  !> Sets the gas that the drops of cell (i,j,k) of block blk see: density rho [kg/m^3], velocity u [m/s], temperature
  !> T [K], gas constant R [J/(kg K)], ratio of specific heats gam, viscosity mu [Pa s] and thermal conductivity kappa
  !> [W/(m K)]. A coupler that owns the gas sets every cell before each ICE step. Writes only that cell.
  subroutine set_gas(blk, i, j, k, rho, u, T, R, gam, mu, kappa)
    use ICE_Advanced_Types_m, only: ICE_block_type
    type(ICE_block_type), intent(inout) :: blk
    integer(I4),          intent(in)    :: i, j, k
    real(R8),             intent(in)    :: rho, u(3), T, R, gam, mu, kappa

    blk%gas_phase%prim(1,i,j,k)   = rho
    blk%gas_phase%prim(2:4,i,j,k) = u
    blk%gas_phase%prim(5,i,j,k)   = T
    blk%gas_phase%R(i,j,k)        = R
    blk%gas_phase%gam(i,j,k)      = gam
    blk%gas_phase%mu(i,j,k)       = mu
    blk%gas_phase%k(i,j,k)        = kappa

  end subroutine set_gas


  !> What family p gives the gas from cell (i,j,k) of block blk, per unit volume and time, over its last step. With
  !> several stages the rates are weighted as the update weighs the stages, so over the step the gas receives exactly
  !> what the family lost. Read it after an ICE step, before the next one.
  pure function exchange_rate(blk, p, i, j, k) result(x)
    use ICE_Advanced_Types_m, only: ICE_block_type
    use ICE_Global_m,         only: nbase
    type(ICE_block_type), intent(in) :: blk
    integer(I4),          intent(in) :: p, i, j, k
    type(exchange_t) :: x
    integer(I4) :: ie

    ! source is what the drops gain; gen_rate is the mass their burning generated, with its momentum and energy
    ie = nbase(p) - 1
    x%generated = blk%cond_phase(p)%gen_rate(1,i,j,k)
    x%mass      = - blk%cond_phase(p)%source(1,i,j,k)   - x%generated
    x%momentum  = - blk%cond_phase(p)%source(2:4,i,j,k) - blk%cond_phase(p)%gen_rate(2:4,i,j,k)
    x%energy    = - blk%cond_phase(p)%source(ie,i,j,k)  - blk%cond_phase(p)%gen_rate(5,i,j,k)
    x%heat      =   blk%cond_phase(p)%heat(i,j,k)

  end function exchange_rate


  !> What the burn-out removed from cell (i,j,k) of family p since the coupler last took it, per unit volume. It is an
  !> amount, not a rate: a time-accurate run removes the last drops of a cell at once, at the end of a step (see
  !> burnout_switch). The cell's record is emptied. Writes only that cell.
  subroutine take_removed(blk, p, i, j, k, x)
    use ICE_Advanced_Types_m, only: ICE_block_type
    use ICE_Global_m,         only: mat_of
    use ICE_Config_Types_m,   only: obj_condensed
    type(ICE_block_type), intent(inout) :: blk
    integer(I4),          intent(in)    :: p, i, j, k
    type(exchange_t),     intent(out)   :: x

    x%mass      = blk%cond_phase(p)%removed(1,i,j,k)
    x%momentum  = blk%cond_phase(p)%removed(2:4,i,j,k)
    x%energy    = blk%cond_phase(p)%removed(5,i,j,k)
    x%generated = blk%cond_phase(p)%removed(6,i,j,k)
    ! all of the removed drops burn
    x%heat      = obj_condensed(mat_of(p))%qComb*(x%mass + x%generated)
    blk%cond_phase(p)%removed(:,i,j,k) = 0._R8

  end subroutine take_removed


  ! ---------------------------------------------------------------------------------------------------------------
  ! Births
  ! ---------------------------------------------------------------------------------------------------------------

  !> Declares that material m generates nu kilograms of the material named `generated` for each kilogram it burns.
  !> The generated mass is born as drops of that material's family (birth_family), with the velocity and temperature
  !> of the burning drops, each drop with the diameter d-gen of m's row in the models file. ICE running alone generates
  !> nothing: a coupler that knows the reaction (from the products file ATLAS writes) declares it. why is empty when
  !> the declaration is taken, else it says why it is not, and nothing is set.
  subroutine set_generation(m, generated, nu, why)
    use ICE_Config_Types_m, only: obj_condensed
    integer(I4),                   intent(in)  :: m
    character(len=*),              intent(in)  :: generated
    real(R8),                      intent(in)  :: nu
    character(len=:), allocatable, intent(out) :: why
    integer(I4) :: p

    why = ''
    if (obj_condensed(m)%combSelect == 0) then
      why = trim(obj_condensed(m)%name)//' does not burn, so it generates nothing'
      return
    endif
    if (.not. nu > 0._R8) then
      why = 'the mass generated per kilogram burnt must be positive'
      return
    endif
    if (.not. obj_condensed(m)%dGen > 0._R8) then
      why = 'the row of '//trim(obj_condensed(m)%name)//' in the models file gives no d-gen, the diameter each '// &
            'generated drop is born with'
      return
    endif
    p = birth_family(generated, why)
    if (p == 0) return
    obj_condensed(m)%genFam = p
    obj_condensed(m)%nuGen  = nu

  end subroutine set_generation


  !> The family in which a mass of the material called `name` is born, or 0 when there is none. why is empty when
  !> there is one, else it says why not. A birth adds mass and a number of drops of one diameter, so the family must
  !> use the MK closure; its material must neither evaporate nor burn; and it must be the material's only family, or
  !> the birth would not know which one to go to.
  function birth_family(name, why) result(p)
    use ICE_Global_m,       only: ngroups, mat_of
    use ICE_Config_Types_m, only: obj_condensed, obj_time_scheme
    character(len=*),              intent(in)  :: name
    character(len=:), allocatable, intent(out) :: why
    integer(I4) :: p
    integer(I4) :: m, q, n

    p   = 0
    why = ''
    m = 0
    do q = 1, size(obj_condensed)
      if (trim(obj_condensed(q)%name) == trim(name)) m = q
    enddo
    if (m == 0) then
      why = trim(name)//' is not a material of ICE''s phase'
      return
    endif
    n = 0
    do q = 1, ngroups
      if (mat_of(q) /= m) cycle
      n = n + 1
      p = q
    enddo
    if (n /= 1) then
      p = 0
      why = trim(name)//' has no family in ICE'
      if (n > 1) why = trim(name)//' has several families in ICE, and a birth goes to one: give it one population'
      return
    endif
    if (trim(obj_time_scheme%model(p)) /= 'MK') then
      why = 'the family of '//trim(name)//' uses the '//trim(obj_time_scheme%model(p))//' closure: a birth adds '// &
            'drops of one diameter, which needs closure = MK'
      p = 0
    elseif (obj_condensed(m)%evapSelect > 0 .or. obj_condensed(m)%combSelect > 0) then
      why = trim(name)//' evaporates or burns, and a mass born in ICE does neither'
      p = 0
    endif

  end function birth_family


  !> The packet of a mass m born in family p with velocity u [m/s] at temperature T [K], every drop with diameter d [m]:
  !> [m, m u, m (e(T) + |u|^2/2), number of drops], with e and the density of the family's material. m is an amount
  !> [kg/m^3] or a rate [kg/(m^3 s)], and the packet is the same.
  pure function birth_packet(p, m, u, T, d) result(pk)
    use ICE_Parameters_m,   only: pi
    use ICE_Global_m,       only: mat_of
    use ICE_Config_Types_m, only: obj_condensed
    use ICE_Lib_Properties, only: mat_rho, mat_e, mat_cp_const
    integer(I4), intent(in) :: p
    real(R8),    intent(in) :: m, u(3), T, d
    real(R8) :: pk(6), e

    if (obj_condensed(mat_of(p))%cs_varies) then
      e = mat_e(obj_condensed(mat_of(p)), T)
    else
      e = mat_cp_const(obj_condensed(mat_of(p)))*T
    endif
    pk(1)   = m
    pk(2:4) = m*u
    pk(5)   = m*(e + 0.5_R8*dot_product(u, u))
    pk(6)   = m/(pi/6._R8*mat_rho(obj_condensed(mat_of(p)), T)*d**3)

  end function birth_packet


  !> A mass born at once in cell (i,j,k) of family p: the packet pk, an amount per unit volume, added to the cell's
  !> conservative state. Writes only that cell.
  subroutine add_birth(blk, p, i, j, k, pk)
    use ICE_Advanced_Types_m, only: ICE_block_type
    use ICE_Global_m,         only: ncond, nbase, mat_of
    use ICE_Config_Types_m,   only: obj_condensed
    use ICE_Lib_Model,        only: model_prim_2_cons, model_cons_2_prim
    type(ICE_block_type), intent(inout) :: blk
    integer(I4),          intent(in)    :: p, i, j, k
    real(R8),             intent(in)    :: pk(6)
    real(R8) :: c(ncond(p))

    c = model_prim_2_cons(p, blk%cond_phase(p)%prim(:,i,j,k), obj_condensed(mat_of(p)))
    c(1:4)        = c(1:4)        + pk(1:4)
    c(nbase(p)-1) = c(nbase(p)-1) + pk(5)
    c(nbase(p))   = c(nbase(p))   + pk(6)
    blk%cond_phase(p)%prim(:,i,j,k) = model_cons_2_prim(p, c, obj_condensed(mat_of(p)))

  end subroutine add_birth


  !> Adds the packet pk, a rate per unit volume and time, to the births in force in cell (i,j,k) of family p. Every
  !> stage of the family's steps receives the births in force over its own step, until clear_birth_rates. Call
  !> clear_birth_rates once first. Writes only that cell.
  subroutine add_birth_rate(blk, p, i, j, k, pk)
    use ICE_Advanced_Types_m, only: ICE_block_type
    type(ICE_block_type), intent(inout) :: blk
    integer(I4),          intent(in)    :: p, i, j, k
    real(R8),             intent(in)    :: pk(6)

    if (.not. allocated(blk%cond_phase(p)%births_ext)) then
      write(*,'(A,I0,A)') ' [ERROR] [ICE::add_birth_rate] family ', p, ': call clear_birth_rates first'
      error stop 1
    endif
    blk%cond_phase(p)%births_ext(:,i,j,k) = blk%cond_phase(p)%births_ext(:,i,j,k) + pk

  end subroutine add_birth_rate


  !> Clears the birth rates in force in every cell of family p on this rank (add_birth_rate). Outside parallel regions.
  subroutine clear_birth_rates(grid, p)
    use ICE_Advanced_Types_m, only: ICE_domain_type
    use ICE_Mod_MPI,          only: is_local_block
    type(ICE_domain_type), intent(inout) :: grid
    integer(I4),           intent(in)    :: p
    integer(I4) :: b

    do b = 1, grid%nb
      if (.not. is_local_block(b)) cycle
      if (.not. allocated(grid%blk(b)%cond_phase(p)%births_ext)) &
        allocate(grid%blk(b)%cond_phase(p)%births_ext(6, grid%blk(b)%dim(1), grid%blk(b)%dim(2), grid%blk(b)%dim(3)))
      grid%blk(b)%cond_phase(p)%births_ext = 0._R8
    enddo

  end subroutine clear_birth_rates

end module ICE_Mod_Exchange
