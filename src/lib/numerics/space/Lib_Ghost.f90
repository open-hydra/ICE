module ICE_Lib_Ghost
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  use ICE_Global_m
  use ICE_Advanced_Types_m
  use ICE_Config_Types_m, only: obj_time_scheme, obj_condensed, condensed_phase_t
  use ICE_Lib_Properties, only : mat_rho
  use ICE_Mod_Metrics, only : delthe
  use ICE_Lib_MK, only : prim_2_cons_MK, cons_2_prim_MK
  use ICE_Lib_IG, only : prim_2_cons_IG, cons_2_prim_IG
  use ICE_Lib_AG, only : prim_2_cons_AG, cons_2_prim_AG, mirror_tensor_AG
  use ICE_Mod_MPI, only : is_local_block
  use ICE_Mod_GhostExchange, only : exchange_ghost_prim

  implicit none
  private
  public :: compute_ghost, fill_second_ghost

contains

  subroutine compute_ghost (grid)
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4) :: i
    integer(kind=I4) :: bm, pm, im, jm, km, fm
    integer(kind=I4) :: ig, jg, kg
    integer(kind=I4) :: bs, is, js, ks, fs
    integer(kind=I4) :: ic, jc, kc
    real(kind=R8)    :: area, normal(1:3), velocity(1:3), veln

    !> Remote cells read below (connection sources, chimera donors) from their owners
    call exchange_ghost_prim(grid)

    !$OMP PARALLEL DEFAULT(NONE), &
    !$OMP SHARED(grid, ncond, mat_of, obj_time_scheme, obj_condensed), &
    !$OMP PRIVATE(i, bm, pm, im, jm, km, fm, ig, jg, kg, bs, is, js, ks, fs, ic, jc, kc, area, normal, velocity, veln)
    !$OMP DO SCHEDULE (dynamic)
    do i = 1, size(grid%bc)

      !> Preliminary assignments
      bm = grid%bc(i)%b
      if (.not. is_local_block(bm)) cycle
      im = grid%bc(i)%i
      jm = grid%bc(i)%j
      km = grid%bc(i)%k
      pm = grid%bc(i)%p
      fm = grid%bc(i)%f

      !> Ghost cell coordinates
      ig = im - guide(fm,1)
      jg = jm - guide(fm,2)
      kg = km - guide(fm,3)

      select case ( grid%bc(i)%type )

        case (101, 201) !> connection / periodic
          bs = grid%bc(i)%bs
          is = grid%bc(i)%is
          js = grid%bc(i)%js
          ks = grid%bc(i)%ks

          grid%blk(bm)%cond_phase(pm)%prim(1:ncond(pm),ig,jg,kg) = grid%blk(bs)%cond_phase(pm)%prim(1:ncond(pm),is,js,ks)


        case (200, 300, 301) !> symmetry (300), the dispersed-phase wall ATLAS tags 301 and the wedge side
                             !>  face (200): a 300/301 face extrapolates particles moving towards it and
                             !>  mirrors the others; 200 and face 3, usually the axis, always mirror
          if (fm <= 2) then
            ic = im - mod(fm,2)
            normal = grid%blk(bm)%dir(1)%f(ic,jm,km)%N
          elseif (fm <= 4) then
            jc = jm - mod(fm,2)
            normal = grid%blk(bm)%dir(2)%f(im,jc,km)%N
          else
            kc = km - mod(fm,2)
            normal = grid%blk(bm)%dir(3)%f(im,jm,kc)%N
          endif

          grid%blk(bm)%cond_phase(pm)%prim(1:ncond(pm),ig,jg,kg) = grid%blk(bm)%cond_phase(pm)%prim(1:ncond(pm),im,jm,km)

          velocity = grid%blk(bm)%cond_phase(pm)%prim(2:4,im,jm,km)
          veln     = dot_product(velocity,normal)
          if (grid%bc(i)%type == 200 .or. veln*real(1-2*mod(fm,2))<=0._R8 .or. fm==3) then
            velocity = velocity - 2._R8*veln*normal
            if (ncond(pm) == 12) grid%blk(bm)%cond_phase(pm)%prim(5:10,ig,jg,kg) = &
              mirror_tensor_AG(grid%blk(bm)%cond_phase(pm)%prim(5:10,im,jm,km), normal)
          endif

          grid%blk(bm)%cond_phase(pm)%prim(2:4,ig,jg,kg) = velocity(1:3)


        case (0) !> planar 2-D face: copy of the interior cell
          grid%blk(bm)%cond_phase(pm)%prim(1:ncond(pm),ig,jg,kg) = grid%blk(bm)%cond_phase(pm)%prim(1:ncond(pm),im,jm,km)


        case (401:403) !> inflow
          if (fm <= 2) then
            ic = im - mod(fm,2)
            normal = grid%blk(bm)%dir(1)%f(ic,jm,km)%N
            area   = grid%blk(bm)%dir(1)%f(ic,jm,km)%A
          elseif (fm <= 4) then
            jc = jm - mod(fm,2)
            normal = grid%blk(bm)%dir(2)%f(im,jc,km)%N
            area   = grid%blk(bm)%dir(2)%f(im,jc,km)%A
          else
            kc = km - mod(fm,2)
            normal = grid%blk(bm)%dir(3)%f(im,jm,kc)%N
            area   = grid%blk(bm)%dir(3)%f(im,jm,kc)%A
          endif

          velocity = grid%blk(bm)%cond_phase(pm)%prim(2:4,im,jm,km)
          veln     = dot_product(velocity,normal)

          !> Switch to extrapolation boundary condition
          if (veln*real(1-2*mod(fm,2))>0._R8) then
            grid%blk(bm)%cond_phase(pm)%prim(1:ncond(pm),ig,jg,kg) = grid%blk(bm)%cond_phase(pm)%prim(1:ncond(pm),im,jm,km)
            cycle
          endif

          !> Switch to symmetry boundary condtion
          if (grid%bc(i)%massflux == 0._R8) then
            grid%blk(bm)%cond_phase(pm)%prim(1:ncond(pm),ig,jg,kg) = grid%blk(bm)%cond_phase(pm)%prim(1:ncond(pm),im,jm,km)
            velocity = velocity - 2._R8*veln*normal
            grid%blk(bm)%cond_phase(pm)%prim(2:4,ig,jg,kg) = velocity(1:3)
            if (ncond(pm) == 12) grid%blk(bm)%cond_phase(pm)%prim(5:10,ig,jg,kg) = &
              mirror_tensor_AG(grid%blk(bm)%cond_phase(pm)%prim(5:10,im,jm,km), normal)


          !> Injection boundary condition
          else

            associate ( prim => grid%blk(bm)%cond_phase(pm)%prim )

            prim(1:ncond(pm),ig,jg,kg) = 0._R8

            if ( ( grid%bc(i)%alpha == 0d0 ) .and. ( grid%bc(i)%beta == 0d0 ) ) then
              grid%bc(i)%alpha = atan2(normal(2),normal(1)+1d-20)
              grid%bc(i)%beta  = atan2(normal(3),normal(1)+1d-20)
              if (abs(normal(3)) < 1d-10) grid%bc(i)%beta = 0.d0
            endif

            !> Assigned gaseous carrier phase
            if (grid%bc(i)%type == 401) then
              veln = dot_product(grid%blk(bm)%gas_phase%prim(2:4,im,jm,km),normal)
              !> Velocity
              prim(2,ig,jg,kg) = grid%bc(i)%velocity * veln * cos(grid%bc(i)%beta) * cos(grid%bc(i)%alpha)
              prim(3,ig,jg,kg) = grid%bc(i)%velocity * veln * cos(grid%bc(i)%beta) * sin(grid%bc(i)%alpha)
              prim(4,ig,jg,kg) = grid%bc(i)%velocity * veln * sin(grid%bc(i)%beta)
              !> Density
              prim(1,ig,jg,kg) = grid%bc(i)%massflux/(1._R8-grid%bc(i)%massflux) * grid%blk(bm)%gas_phase%prim(1,im,jm,km) / grid%bc(i)%velocity
              !> Temperature
              prim(ncond(pm)-1,ig,jg,kg) = grid%bc(i)%temperature * grid%blk(bm)%gas_phase%prim(5,im,jm,km)

            !> Direct assignement of massflux, velocity, and temperature.
            elseif (grid%bc(i)%type == 402) then
              !> Velocity
              prim(2,ig,jg,kg) = grid%bc(i)%velocity * cos(grid%bc(i)%beta) * cos(grid%bc(i)%alpha)
              prim(3,ig,jg,kg) = grid%bc(i)%velocity * cos(grid%bc(i)%beta) * sin(grid%bc(i)%alpha)
              prim(4,ig,jg,kg) = grid%bc(i)%velocity * sin(grid%bc(i)%beta)
              !> Density
              veln = dot_product(prim(2:4,ig,jg,kg),normal)
              prim(1,ig,jg,kg) = grid%bc(i)%massflux / abs(veln)
              !> Temperature
              prim(ncond(pm)-1,ig,jg,kg) = grid%bc(i)%temperature

            !> Direct assignement of massflux and temperature. Velocity is computed from the gas phase
            elseif (grid%bc(i)%type == 403) then
              !> Velocity
              veln = norm2(grid%blk(bm)%gas_phase%prim(2:4,im,jm,km))
              prim(2,ig,jg,kg) = grid%bc(i)%velocity * veln * cos(grid%bc(i)%beta) * cos(grid%bc(i)%alpha)
              prim(3,ig,jg,kg) = grid%bc(i)%velocity * veln * cos(grid%bc(i)%beta) * sin(grid%bc(i)%alpha)
              prim(4,ig,jg,kg) = grid%bc(i)%velocity * veln * sin(grid%bc(i)%beta)
              !> Density
              veln = dot_product(prim(2:4,ig,jg,kg),normal)
              prim(1,ig,jg,kg) = grid%bc(i)%massflux / abs(veln)
              !> Temperature
              prim(ncond(pm)-1,ig,jg,kg) = grid%bc(i)%temperature

            endif

            !> N particles, with the condensed density at the inlet temperature
            prim(ncond(pm),ig,jg,kg) =  prim(1,ig,jg,kg) / mat_rho(obj_condensed(mat_of(pm)), prim(ncond(pm)-1,ig,jg,kg)) / &
                                        (4._R8/3._R8*pi*grid%bc(i)%radius**3._I4)

            !> Pseudo pressure
            select case (trim(obj_time_scheme%model(pm)))
            case ('IG')
              prim(5,ig,jg,kg) = 1d-6
            case ('AG')
              prim(5,ig,jg,kg) = 1d-6
              prim(8,ig,jg,kg) = 1d-6
              prim(10,ig,jg,kg) = 1d-6
            end select

          endassociate

          endif


        case (400) !> extrapolation
          grid%blk(bm)%cond_phase(pm)%prim(1:ncond(pm),ig,jg,kg) = grid%blk(bm)%cond_phase(pm)%prim(1:ncond(pm),im,jm,km)


        case (102) !> chimera: first ghost layer from its donors
          call ghost_chimera(grid%blk, grid%bc(i), 1, obj_condensed(mat_of(pm)))


      end select

    enddo
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine compute_ghost


  subroutine fill_second_ghost(grid)
    implicit none
    type(ICE_domain_type), intent(inout) :: grid
    integer(kind=I4) :: i
    integer(kind=I4) :: bm, pm, im, jm, km, fm
    integer(kind=I4) :: ig, jg, kg
    integer(kind=I4) :: ig2, jg2, kg2
    integer(kind=I4) :: ip, jp, kp
    integer(kind=I4) :: bs, is, js, ks, fs
    integer(kind=I4) :: ic, jc, kc
    real(kind=R8)    :: normal(1:3), velocity(1:3)

    !$OMP PARALLEL DEFAULT(NONE), &
    !$OMP SHARED(grid, ncond, mat_of, obj_condensed), &
    !$OMP PRIVATE(i, bm, pm, im, jm, km, fm, ig, jg, kg, ig2, jg2, kg2, ip, jp, kp, bs, is, js, ks, fs, &
    !$OMP         ic, jc, kc, normal, velocity)
    !$OMP DO SCHEDULE(dynamic)
    do i = 1, size(grid%bc)

      if (grid%bc(i)%type == 0) cycle

      bm = grid%bc(i)%b
      if (.not. is_local_block(bm)) cycle
      im = grid%bc(i)%i ; jm = grid%bc(i)%j ; km = grid%bc(i)%k
      pm = grid%bc(i)%p ; fm = grid%bc(i)%f

      ig  = im -   guide(fm,1) ; jg  = jm -   guide(fm,2) ; kg  = km -   guide(fm,3)
      ig2 = im - 2*guide(fm,1) ; jg2 = jm - 2*guide(fm,2) ; kg2 = km - 2*guide(fm,3)

      select case (grid%bc(i)%type)

        case (101, 201) !> connection: copy second interior cell of source block
          bs = grid%bc(i)%bs ; fs = grid%bc(i)%fs
          is = grid%bc(i)%is + guide(fs,1)
          js = grid%bc(i)%js + guide(fs,2)
          ks = grid%bc(i)%ks + guide(fs,3)
          grid%blk(bm)%cond_phase(pm)%prim(1:ncond(pm),ig2,jg2,kg2) = &
            grid%blk(bs)%cond_phase(pm)%prim(1:ncond(pm),is,js,ks)

        case (102) !> chimera: second ghost layer from its own donors
          call ghost_chimera(grid%blk, grid%bc(i), 2, obj_condensed(mat_of(pm)))

        case (200) !> wedge side face: the second ghost mirrors the second interior cell
          ip = im + guide(fm,1) ; jp = jm + guide(fm,2) ; kp = km + guide(fm,3)
          if (fm <= 2) then
            ic = im - mod(fm,2)
            normal = grid%blk(bm)%dir(1)%f(ic,jm,km)%N
          elseif (fm <= 4) then
            jc = jm - mod(fm,2)
            normal = grid%blk(bm)%dir(2)%f(im,jc,km)%N
          else
            kc = km - mod(fm,2)
            normal = grid%blk(bm)%dir(3)%f(im,jm,kc)%N
          endif
          grid%blk(bm)%cond_phase(pm)%prim(1:ncond(pm),ig2,jg2,kg2) = grid%blk(bm)%cond_phase(pm)%prim(1:ncond(pm),ip,jp,kp)
          velocity = grid%blk(bm)%cond_phase(pm)%prim(2:4,ip,jp,kp)
          grid%blk(bm)%cond_phase(pm)%prim(2:4,ig2,jg2,kg2) = velocity - 2._R8*dot_product(velocity,normal)*normal
          if (ncond(pm) == 12) grid%blk(bm)%cond_phase(pm)%prim(5:10,ig2,jg2,kg2) = &
            mirror_tensor_AG(grid%blk(bm)%cond_phase(pm)%prim(5:10,ip,jp,kp), normal)

        case default !> 2nd-order extrapolation: P(g2) = 3*P(g1) - 3*P(m) + P(m+1)
          ip = im + guide(fm,1) ; jp = jm + guide(fm,2) ; kp = km + guide(fm,3)
          grid%blk(bm)%cond_phase(pm)%prim(1:ncond(pm),ig2,jg2,kg2) = &
            3d0 * grid%blk(bm)%cond_phase(pm)%prim(1:ncond(pm),ig,jg,kg) - &
            3d0 * grid%blk(bm)%cond_phase(pm)%prim(1:ncond(pm),im,jm,km) + &
                  grid%blk(bm)%cond_phase(pm)%prim(1:ncond(pm),ip,jp,kp)

      end select

    enddo
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine fill_second_ghost


  !> Chimera ghost cell of layer g: conservative blend of the donor cells found by
  !> ATLAS BCB, cons(ghost) = sum_c w_c * cons(donor_c), as MOSE's Ghost_Chimera does.
  !> Donors are interior cells, so this only reads interior data and is safe inside
  !> the OMP loops over bc entries.
  subroutine ghost_chimera(blk, bc, g, mat)
    implicit none
    type(ICE_block_type), intent(inout) :: blk(:)
    type(ICE_bc_type),    intent(in)    :: bc
    integer(kind=I4),     intent(in)    :: g
    type(condensed_phase_t), intent(in) :: mat
    integer(kind=I4) :: c, c1, c2, pm, nv, bs, is, js, ks, ig, jg, kg
    real(kind=R8)    :: consg(12)

    pm = bc%p
    nv = ncond(pm)
    if (g == 1) then
      c1 = 1          ; c2 = bc%ni(1)
    else
      c1 = bc%ni(1)+1 ; c2 = sum(bc%ni)
    endif

    consg = 0._R8
    do c = c1, c2
      bs = bc%donorID(c,1) ; is = bc%donorID(c,2) ; js = bc%donorID(c,3) ; ks = bc%donorID(c,4)
      consg(1:nv) = consg(1:nv) + bc%volume_fraction(c) * &
                    model_prim_2_cons(pm, blk(bs)%cond_phase(pm)%prim(1:nv,is,js,ks), mat)
    enddo

    ig = bc%i - g*guide(bc%f,1) ; jg = bc%j - g*guide(bc%f,2) ; kg = bc%k - g*guide(bc%f,3)
    blk(bc%b)%cond_phase(pm)%prim(1:nv,ig,jg,kg) = model_cons_2_prim(pm, consg(1:nv), mat)

  end subroutine ghost_chimera


  !> Group-specific conversions. The model procedure pointers in ICE_Lib_Model are bound to
  !> one group at a time, while the ghost loops cover every group, so dispatch explicitly.
  function model_prim_2_cons(p, prim, mat) result(cons)
    implicit none
    integer(kind=I4), intent(in) :: p
    real(kind=R8),    intent(in) :: prim(:)
    type(condensed_phase_t), intent(in) :: mat
    real(kind=R8)                :: cons(size(prim))

    select case (trim(obj_time_scheme%model(p)))
    case ('MK'); cons = prim_2_cons_MK(prim, mat)
    case ('IG'); cons = prim_2_cons_IG(prim, mat)
    case ('AG'); cons = prim_2_cons_AG(prim, mat)
    end select

  end function model_prim_2_cons


  function model_cons_2_prim(p, cons, mat) result(prim)
    implicit none
    integer(kind=I4), intent(in) :: p
    real(kind=R8),    intent(in) :: cons(:)
    type(condensed_phase_t), intent(in) :: mat
    real(kind=R8)                :: prim(size(cons))

    select case (trim(obj_time_scheme%model(p)))
    case ('MK'); prim = cons_2_prim_MK(cons, mat)
    case ('IG'); prim = cons_2_prim_IG(cons, mat)
    case ('AG'); prim = cons_2_prim_AG(cons, mat)
    end select

  end function model_cons_2_prim


end module ICE_Lib_Ghost
