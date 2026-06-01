module ICE_Lib_Ghost
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  use ICE_Global_m
  use ICE_Advanced_Types_m
  use ICE_Config_Types_m, only: obj_condensed, obj_time_scheme
  use ICE_Mod_Metrics, only : delthe

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

    !$OMP PARALLEL DEFAULT(NONE), &
    !$OMP SHARED(grid, ncond, obj_condensed, obj_time_scheme), &
    !$OMP PRIVATE(i, bm, pm, im, jm, km, fm, ig, jg, kg, bs, is, js, ks, fs, ic, jc, kc, area, normal, velocity, veln)
    !$OMP DO SCHEDULE (dynamic)
    do i = 1, size(grid%bc)

      !> Preliminary assignments
      bm = grid%bc(i)%b
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

        case (1) !> connection
          bs = grid%bc(i)%bs
          is = grid%bc(i)%is
          js = grid%bc(i)%js
          ks = grid%bc(i)%ks

          grid%blk(bm)%cond_phase(pm)%prim(1:ncond(pm),ig,jg,kg) = grid%blk(bs)%cond_phase(pm)%prim(1:ncond(pm),is,js,ks)


        case (3, 5, 6) !> wall - Use extrapolation when the particles are moving towards the wall. Otherwise, symmetry
                       !>        Symmetry is enforced on face 3 which is usually the symmetry axis
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
          if (veln*real(1-2*mod(fm,2))<=0._R8 .or. fm==3) then
            velocity = velocity - 2._R8*veln*normal
          endif

          grid%blk(bm)%cond_phase(pm)%prim(2:4,ig,jg,kg) = velocity(1:3)


        case (4,14) !> inflow
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
            if (grid%bc(i)%injtype == 0) then
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
            elseif (grid%bc(i)%injtype == 1) then
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
            elseif (grid%bc(i)%injtype == 2) then
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

            !> N particles
            prim(ncond(pm),ig,jg,kg) =  prim(1,ig,jg,kg) / obj_condensed%rho_al / (4._R8/3._R8*pi*grid%bc(i)%radius**3._I4)

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


        case (11) !> extrapolation
          grid%blk(bm)%cond_phase(pm)%prim(1:ncond(pm),ig,jg,kg) = grid%blk(bm)%cond_phase(pm)%prim(1:ncond(pm),im,jm,km)


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

    !$OMP PARALLEL DEFAULT(NONE), &
    !$OMP SHARED(grid, ncond), &
    !$OMP PRIVATE(i, bm, pm, im, jm, km, fm, ig, jg, kg, ig2, jg2, kg2, ip, jp, kp, bs, is, js, ks, fs)
    !$OMP DO SCHEDULE(dynamic)
    do i = 1, size(grid%bc)

      if (grid%bc(i)%type == 0 .or. grid%bc(i)%type == 2) cycle

      bm = grid%bc(i)%b
      im = grid%bc(i)%i ; jm = grid%bc(i)%j ; km = grid%bc(i)%k
      pm = grid%bc(i)%p ; fm = grid%bc(i)%f

      ig  = im -   guide(fm,1) ; jg  = jm -   guide(fm,2) ; kg  = km -   guide(fm,3)
      ig2 = im - 2*guide(fm,1) ; jg2 = jm - 2*guide(fm,2) ; kg2 = km - 2*guide(fm,3)

      select case (grid%bc(i)%type)

        case (1) !> connection: copy second interior cell of source block
          bs = grid%bc(i)%bs ; fs = grid%bc(i)%fs
          is = grid%bc(i)%is + guide(fs,1)
          js = grid%bc(i)%js + guide(fs,2)
          ks = grid%bc(i)%ks + guide(fs,3)
          grid%blk(bm)%cond_phase(pm)%prim(1:ncond(pm),ig2,jg2,kg2) = &
            grid%blk(bs)%cond_phase(pm)%prim(1:ncond(pm),is,js,ks)

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


end module ICE_Lib_Ghost
