module driver
!use parameters
implicit none

contains

! Numerical functions

!!!! F2PY does not allow for allocatable function inputs!!!
!!!! major limitation..... !!!
subroutine convolve(nb, nx, ny, bb, xx, yy)
    integer, intent(in) :: nb             !# number of coefficients in filter
    integer, intent(in) :: nx             !# number of coefficients in input
    integer, intent(in) :: ny             !# number of coefficients in output will be nx+nb-1
    real, intent(in), dimension(nb) :: bb         !# filter coefficients
    real, intent(in), dimension(nx) :: xx         !# input trace
    real, intent(out), dimension(ny) :: yy

    !internal
    integer :: ib, ix, iy
    !integer, dimension(ny) :: yy
    yy = 0.0
    !ny = nx + nb -1
    do ib = 1, nb
        do ix = 1, nx
            yy(ix+ib-1) = yy( ix+ib-1) + xx(ix) * bb(ib)
        end do
    end do
end subroutine convolve


! This is the gamma function kernel used to convolve
! runoff signal to produce streamfow
! WHERE DOES THIS FUNCTION COME FROM?
real function ht(t, k, n) !=3.5, N=4)
    real, intent(in) :: t, k, N
    ht = (t/k)**(N-1) * EXP(-t/k)/k*GAMMA(N)
end function ht

! Helper functions
function PET(L, t)
    real :: PET, L, t
    real :: sat_vap_pres, abs_humidity
    ! Potential evapotranspiration calculator
    ! L: Length of day
    ! rho_v_sat: staturated absolute humidity
    ! etpar: parameter
    sat_vap_pres = (.6112)*exp((17.67*t)/(t + 243.5))
    abs_humidity = (2165 * sat_vap_pres)/(t + 273.15)
    PET = .55 * (L/12.)**2 * abs_humidity
    ! # Potential evapotranspiration calculator
    ! # L: Length of day
    ! # rho_v_sat: staturated absolute humidity
    ! # etpar: parameter
    ! sat_vap_pres = (.6112)*np.exp((17.67*t)/(t + 243.5))
    ! abs_humidity = (2165. * sat_vap_pres)/(t + 273.15)

    ! # compute PET
    ! PET=.55 * (L/12.)**2 * abs_humidity
    ! return PET
end function PET



subroutine model_driver(HYDROOP,    &         ! OPTION Hydrologic Model Option
                        ntimes,     &         ! FORCING   Number of model timesteps
                        PET,        &         ! FORCING   Potential Evapotranspiration
                        tair,       &         ! FORCING   Air Temperature, length N
                        qin,        &         ! FORCING   Air Temperature, length N
                        sm1i,       &         ! INTITIAL SM CONDITION
                        sm2i,       &         ! INTITIAL SM CONDITION
                        sm1max,     &         ! PARAMETER   ?
                        sm2max,     &         ! PARAMETER   ?
                        fracten,   &         ! PARAMETER   ?
                        ku,         &         ! PARAMETER   PERCOLATION
                        c,          &         ! PARAMETER   PERCOLATION
                        alpha,      &         ! PARAMETER   PERCOLATION   --- OPTIONAL
                        ks,         &         ! PARAMETER   BASEFLOW      --- OPTIONAL
                        lam,        &         ! PARAMETER   BASEFLOW      --- OPTIONAL
                        lowercasen, &         ! PARAMETER   BASEFLOW      --- OPTIONAL
                        beta,       &         ! PARAMETER   SFROFF        --- OPTIONAL
                        Nr,         &         !
                        kr,         &         !
                        chanVecOutput, eVecOutput, smOutput)          ! OUTPUT


    !use snowmodule17

    ! MODEL RUN OPTIONS
    integer, intent(in) :: HYDROOP


    ! MODEL FORCINGS
    integer, intent(in) :: ntimes
    real, intent(in), dimension(ntimes) :: PET
    real, intent(in), dimension(ntimes) :: tair
    real, intent(in), dimension(ntimes) :: qin ! change this later. for testing purposes

    !SOIL PARAMETERS
    real, intent(in) :: sm1i
    real, intent(in) :: sm2i
    real, intent(in) :: sm1max
    real, intent(in) :: sm2max

    ! MORE SOIL PARAMETERS
    real, intent(in) :: fracten

    ! PERCOLATION PARAMETERS
    real, intent(in), optional :: ku       ! Percolation option A,?
    real, intent(in), optional :: c        ! Percolation option A,? exponent
    real, intent(in), optional :: alpha    ! Percolation option C

    ! BASEFLOW PARAMETERS
    real, intent(in), optional :: ks
    real, intent(in), optional :: lam
    real, intent(in), optional :: lowercasen

    ! SFROFF PARAMETERS
    real, intent(in), optional :: beta     ! Option 0 -- ?



    ! ROUTING PARAMETERS
    real, intent(in) :: Nr
    real, intent(in) :: kr

    ! OUTPUT
    real, intent(out), dimension(ntimes) :: chanVecOutput  ! streamflow after routing
    real, intent(out), dimension(ntimes) :: eVecOutput
    real, intent(out), dimension(ntimes) :: smOutput

    !real, intent(out), dimension(nlayers,ntimes) :: sweVecOutput

    ! INTERNAL
    integer :: i                     ! timestep
    integer :: ny
    real                        :: sm1tmax          ! max tension storage top layer
    real                        :: sm2tmax          ! max tension storage bottom layer
    real, dimension(ntimes)     :: qVecOutput       ! output streamflow
    real, dimension(ntimes)     :: sm1              ! state (sm in top layer)
    real, dimension(ntimes)     :: sm2              ! state (sm in bottom layer)
    real, dimension(ntimes)     :: sm1t             ! state (free water in sm1?)
    real, dimension(ntimes)     :: sm2t             ! state (free water in sm1?)
    real, dimension(ntimes)     :: e1               ! flux from sm1 --> atmos. (evaporation)
    real, dimension(ntimes)     :: e2               ! flux from sm1 --> atmos. (evaporation)
    real, dimension(ntimes)     :: q12              ! flux from sm1 --> sm2
    real, dimension(ntimes)     :: qdir             ! flux from sm1 --> sm2
    real, dimension(ntimes)     :: qb               ! flux from sm2 --> channel (baseflow )
    real, dimension(ntimes)     :: qsx              ! Overland flow
    real, dimension(ntimes)     :: overflow_sm1     ! Bucket overflow
    real, dimension(ntimes)     :: overflow_sm2     ! Bucket overflow
    real, dimension(ntimes)     :: deficit_sm1      ! Bucket overflow
    real, dimension(ntimes)     :: deficit_sm2      ! Bucket overflow
    real, dimension(ntimes)     :: htv              ! Routing kernel
    real, dimension(ntimes*2-1) :: cnvrt

    ! Scalars...
    real :: ac                      ! Saturated area
    real :: pet_remainder           ! unsatisfied PET from the 1st layer

! Collapse Melt into the single step
!    real, dimension(ntimes) :: outflowVecTotal
!    real, dimension(ntimes) :: outflowVec_adjusted
!    real :: scale ! corrects runoff

    real :: residual, ds1, ds2, hold1, hold2
    real, parameter:: tol=.1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!           INITIAL CONDITIONS/BNDARDY                   !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    qVecOutput = 0.
    eVecOutput = 0.
    chanVecOutput = 0.
    smOutput = 0.


    qsx = 0.
    e1 = 0.
    e2 = 0.
    sm1 = 0.
    sm2 = 0.
    qdir = 0.
    qb  = 0.
    q12 = 0.
    overflow_sm1 = 0.
    overflow_sm2 = 0.
    pet_remainder = 0.

    ! initial conditions...
    sm1(1) = sm1i*sm1max
    sm2(1) = sm2i*sm2max


    ! qin() = 0.
    ! PET(ntimes) = 0.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!           BEGIN MAIN MODEL EXECUTION           !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Begin time stepping (1st order explicit Runge-Kutta)
    ! ET, soilwater, flow must be called simultaneously
    ! since they all depend on soil water values
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    select case(HYDROOP)

        ! THESE ARE THE DIFFERENT HYDRO MODEL OPTIONS



        ! BEGIN FIRST MODEL OPTION -- ARNO/VIC-like model...
        case(0)

        sm1tmax = sm1max * fracten
        sm2tmax = sm2max * fracten
        do i=1,ntimes-1

            ! Compute fraction of water that is stored in tension
            sm1t(i) = MIN(sm1(i), sm1tmax)

            ! Compute ET
            !e1(i) = PET(i) * (sm1(i)/sm1max)**etpar
            e1(i) = PET(i) * MIN(sm1t(i), sm1tmax)/sm1tmax

            if (e1(i) > sm1(i)) then
                  e1(i) = sm1(i)
            end if

            ! ! subtract it from sm1
            sm1(i) = sm1(i) - e1(i)

            ! Compute saturated area
            ac = MAX(0., 1. - (1. -  sm1(i)/sm1max) ** beta)

            ! Compute surface runoff
            qsx(i) = ac*qin(i)

            ! Compute the "direct" precipitation component
            qdir(i) = qin(i) - qsx(i)

            ! Add qdir to the top bucket
            sm1(i) = sm1(i) + qdir(i)

            ! remove the overflow
            if (sm1(i) > sm1max) then
                overflow_sm1(i)  = sm1(i) - sm1max
                sm1(i) = sm1(i) - overflow_sm1(i)
            else
                overflow_sm1(i) = 0.0
            end if

            ! Now we add the percolation to the lower layer ... compute it first
            q12(i) = ku * (sm1(i)/sm1max)**c
            if (q12(i) > sm1(i)) then
                  q12(i) = sm1(i)
            end if
            sm1(i) = sm1(i) - q12(i)

            ! ! ! -- Now do lower layer stuff ----

            ! Compute baseflow and remove it
            qb(i) = ks*(sm2(i)/sm2max)**lowercasen
            if (qb(i) > sm2(i)) then
                  qb(i) = sm2(i)
            end if
            sm2(i) = sm2(i) - qb(i)

            ! add in the percolation component
            sm2(i) = sm2(i) + q12(i)

            if (sm2(i) > sm2max) then
                overflow_sm2(i) = sm2(i) - sm2max
                sm2(i) = sm2(i) - overflow_sm2(i)
            else
                overflow_sm2(i) = 0.0
            end if


            ! calculate the residual lower layer ET first
            pet_remainder =  max(0., PET(i) - e1(i))

            ! compute how much water is in tension storage
            sm2t(i) = min(sm2(i), sm2tmax)

            ! Now compute the residual from the lower bucket
            !e2(i) = pet_remainder*(sm2(i)/sm2max)
            e2(i) = pet_remainder*min(sm2t(i), sm2tmax)/sm2tmax

            if (e2(i) > sm2(i)) then
                  e2(i) = sm2(i)
            end if

            ! ! ! remove the ET from the lower layer
            sm2(i) = sm2(i) - e2(i)

            ! ! ! update the next states...
            sm1(i+1) = sm1(i)
            sm2(i+1) = sm2(i)

            ! ! store streamflow
            qVecOutput(i) = qb(i) + qsx(i) + overflow_sm1(i) + overflow_sm2(i)
            eVecOutput(i) = e1(i) + e2(i)
            smOutput(i) = sm1(i) + sm2(i)
        end do
        ! END FIRST MODEL OPTION

        ! -----------
        ! BEGIN SECOND MODEL OPTION
        ! case(1)
        ! do i=1,ntimes-1

        !     !!!!!!! ET OPTION --- Eagelson !!!!!
        !     ! sm1_crit is the sm where veg closes stomata
        !     ! sm1_wilt is the sm where plants begin to wilt
        !     ! fveg = max( min((sm1 - sm1_wilt)/(sm1_crit - sm1_wilt), 1), 0)

        !     ! ! Compute ET
        !     ! e1(i) = PET(i) * (sm1(i)/sm1max)

        !     ! !




        !     ! ! subtract it from sm1
        !     sm1(i) = sm1(i) - e1(i)

        !     ! Compute saturated area
        !     ac = max(0., 1. - (1. -  sm1(i)/sm1max) ** beta)

        !     ! Compute surface runoff
        !     qsx(i) = ac*qin(i)

        !     ! Compute the "direct" precipitation component
        !     qdir(i) = qin(i) - qsx(i)

        !     ! Add qdir to the top bucket
        !     sm1(i) = sm1(i) + qdir(i)

        !     ! remove the overflow
        !     if (sm1(i) > sm1max) then
        !         overflow_sm1(i)  = sm1(i) - sm1max
        !         sm1(i) = sm1(i) - overflow_sm1(i)
        !     else
        !         overflow_sm1(i) = 0.0
        !     end if

        !     ! Now we add the percolation to the lower layer ... compute it first
        !     q12(i) = ku * (sm1(i)/sm1max)**c
        !     if (q12(i) > sm1(i)) then
        !           q12(i) = sm1(i)
        !     end if
        !     sm1(i) = sm1(i) - q12(i)

        !     ! ! ! -- Now do lower layer stuff ----

        !     ! Compute baseflow and remove it
        !     qb(i) = ks*(sm2(i)/sm2max)**lowercasen
        !     if (qb(i) > sm2(i)) then
        !           qb(i) = sm2(i)
        !     end if
        !     sm2(i) = sm2(i) - qb(i)

        !     ! add in the percolation component
        !     sm2(i) = sm2(i) + q12(i)

        !     if (sm2(i) > sm2max) then
        !         overflow_sm2(i) = sm2(i) - sm2max
        !         sm2(i) = sm2(i) - overflow_sm2(i)
        !     else
        !         overflow_sm2(i) = 0.0
        !     end if

        !     ! calculate the residual lower layer ET first
        !     pet_remainder =  max(0., PET(i) - e1(i))

        !      ! Now compute the residual from the lower bucket
        !     e2(i) = pet_remainder*(sm2(i)/sm2max)

        !     if (e2(i) > sm2(i)) then
        !           e2(i) = sm2(i)
        !     end if

        !     ! ! ! remove the ET from the lower layer
        !     sm2(i) = sm2(i) - e2(i)

        !     ! ! ! update the next states...
        !     sm1(i+1) = sm1(i)
        !     sm2(i+1) = sm2(i)

        !     ! ! store streamflow
        !     qVecOutput(i) = qb(i) + qsx(i) + overflow_sm1(i) + overflow_sm2(i)
        !     eVecOutput(i) = e1(i) + e2(i)
        !     smOutput(i) = sm1(i) + sm2(i)
        ! end do


    end select

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!! BEGIN ROUTING MODULE !!!!!!!!!!!!!!!!!!!!!
    !!!     This is where the runoff from the bucket model        !!!
    !!!     gets routed using a unit hydrograph approach.         !!!
    !!!     this means that a gamma function gets convolved       !!!
    !!!     with the spiky-looking runoff timeseries to smooth it !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! ROUTING MODULE
    ! now convert qVecOutput to streamflow...
    do i=1,ntimes
        htv(i) = ht(REAL(i),Nr,Kr)
    end do

    ! normalize it by 1 -- so as to not add/subtract flow
    htv = htv/SUM(htv)
    ny = ntimes*2-1
    call convolve(ntimes, ntimes, ny, htv, qVecOutput, cnvrt)

    ! the convolution adds on a bunch of points at hte end that we don't need
    do i=1,ntimes
        chanVecOutput(i) = cnvrt(i)
    end do

    ! DO SOME CHECKING
    ds1 = sm1(ntimes) - sm1i
    ds2 = sm2(ntimes) - sm2i
    hold1 = SUM(qsx) + SUM(overflow_sm1) + SUM(q12) + SUM(e1)
    hold2 = SUM(qb) + SUM(overflow_sm2) + sum(e2)

    residual = (sum(qin) - ds1 - hold1) + (sum(q12) - ds2 - hold2)
    ! if (ABS(residual) > tol) then
    !      print*, "Water balance error; residual = " , residual
    ! end if

!  print*, params%nmf

end subroutine model_driver


end module driver





