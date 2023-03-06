module snowmodule17
implicit none
contains

!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1) Snow-17 accumulation and ablation model.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Adapted from
! This version of Snow-17 is intended for use at a point location.
! Based on Anderson (2006) and Mark Raleigh's matlab code.
! Primary Citations:
! 1.  Anderson, E. A. (1973), National Weather Service River Forecast System
!     Snow   Accumulation   and   Ablation   Model,   NOAA   Tech.   Memo.   NWS
!     HYDro-17, 217 pp., U.S. Dep. of Commer., Silver Spring, Md.
! 2.  Anderson, E. A. (1976), A point energy and mass balance model of a snow
!     cover, NOAA Tech. Rep. 19, 150 pp., U.S. Dep. of Commer., Silver Spring, Md.
!
! Python code Written by Joe Hamman April, 2013
!
! Translated into Fortran (lol why) by William Rudisill, March 2021
!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real function melt_function(jday,  &
	                          dt,    &
	                          mfmax, &
	                          mfmin)
    implicit none

    ! input
    real, intent(in) :: jday
    real, intent(in) :: dt
    real, intent(in) :: mfmax
    real, intent(in) :: mfmin

    ! internal
    integer :: n_mar21
    integer :: days
    real :: sv
    real :: av
    ! BEGIN
    ! ------

    ! a few parameters
    days = 365.
    n_mar21 = jday - 80.
    sv = (0.5 * SIN((n_mar21 * 2. * 3.14159265359) / days)) + 0.5

    if ((jday <= 77).or.(jday >=227)) then
        ! assumes lat < 45
        av = 0.0

    else if ((jday >= 117).and.(jday <= 227)) then
        av = 1.

    else if ((jday >= 78).and.(jday <=116)) then
        ! interpolate between them... 0 -- 1
        av = (jday-78)/(116.-78.)

    else if ((jday >= 228).and.(jday <= 266)) then
        ! interpolate between them... 0 -- 1
        av = (jday-228)/(266.-228.)

    else !there shouldnt be another else...
        av = 1.
    end if
    ! do something different for northern latitudes... NOT IMPLEMENTED
    ! av = <code here>
    ! ....
    melt_function = (dt / 6) * ((sv * av * (mfmax - mfmin)) + mfmin)
end function melt_function



! Begin subroutines
subroutine snow17(jday, &
  	              precip, &
  	              tair, &
  	              elevation, &
  	              dt, &
  	              rvs, &
  	              uadj, &
  	              mbase, &
  	              mfmax, &
  	              mfmin, &
  	              tipm, &
  	              nmf, &
                  plwhc, &
                  pxtemp, &
                  pxtemp1, &
                  pxtemp2, &
                  swe, &
                  ait, &
                  w_q, &
                  w_i, &
                  deficit, &
                  outflow, &  !ouput
                  rain)       !output

    implicit none

    ! input
    real, intent(in) :: jday
    real, intent(in) :: precip
    real, intent(in) :: tair
    real, intent(in) :: elevation
    real, intent(in) :: dt
    integer, intent(in) :: rvs
    real, intent(in) :: uadj
    real, intent(in) :: mbase
    real, intent(in) :: mfmax
    real, intent(in) :: mfmin
    real, intent(in) :: tipm
    real, intent(in) :: nmf
    real, intent(in) :: plwhc
    real, intent(in) :: pxtemp
    real, intent(in) :: pxtemp1
    real, intent(in) :: pxtemp2

    ! inouts (snow states. these can change)
    real, intent(inout) :: swe
    real, intent(inout) :: ait
    real, intent(inout) :: w_q
    real, intent(inout) :: w_i
    real, intent(inout) :: deficit

    ! OUTPUT VARIABLES-- snowmelt and rain
    real, intent(out) :: outflow
    real, intent(out) :: rain


    !!!!! Declare f2py outputs... inout will not return anything in python
    !!!!! ----------------------------------------------------------------
    !!!! f2py intent(in,out) :: swe    (only implemented w/ one comment )
    !!!! f2py intent(in,out) :: ait    (only implemented w/ one comment )
    !!!! f2py intent(in,out) :: w_qx    (only implemented w/ one comment )
    !!!! f2py intent(in,out) :: w_q    (only implemented w/ one comment )
    !!!! f2py intent(in,out) :: w_i    (only implemented w/ one comment )
    !!!! f2py intent(in,out) :: deficit    (only implemented w/ one comment )

    ! internal
    real, dimension(2) :: transitionx
    real, dimension(2) :: transitiony
    real :: stefan        ! constants
    real :: p_atm
    real :: mf
    real :: tipm_dt
    real :: pn
    real :: w_qx
    real :: delta_hd_snow
    real :: fracrain
    real :: fracsnow
    real :: t_snow_new
    real :: t_rain
    real :: e_sat
    real :: delta_hd_t
    real :: m_ros
    real :: m_ros1
    real :: m_ros2
    real :: m_ros3
    real :: m_nr
    real :: qw    ! liquid water held by snow
    real :: melt  ! snow melt
    real :: ex_lw ! excess liqud water
    real :: k

    !!!!!    !f2py intent(in,out) :: e_sat

    ! BEGIN
    ! ------

    ! Deal with parameters that should not be zero...
    ! if (uadj < 0) then
    !   call exit(0)
    ! end if


    ! constant
    !stefan = 6.12 * (10 ** (-10))
    stefan=6.12E-10

    ! atmospheric pressure (mb) where elevation is in HUNDREDS of meters
    ! (this is incorrectly stated in the manual)
    p_atm = 33.86 * (29.9 - (0.335 * elevation / 100) + (0.00022 * ((elevation / 100) ** 2.4)))

    ! temperature index constant used to compute snow temperature index (ait).
    ! notice this only depends on the input parameters tipm and dt
    tipm_dt = 1.0 - ((1.0 - tipm) ** (dt / 6))

    ! compute melt function
    mf = melt_function(jday, dt, mfmax, mfmin)

    ! saturated vapor pressure at tair (mb)
    e_sat = 2.7489 * (10 ** 8) * EXP((-4278.63 / (tair + 242.792)))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Select the 'rain versus snow' method !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Select rain vs snow method
    select case(rvs)
      ! Option 0
      case(0);
        if (tair <= pxtemp) then
            ! then the air temperature is cold enough for snow to occur
            fracsnow = 1.0
        else
            fracsnow = 0.0
        end if

      ! Option 1
      case(1)
         if (tair <= pxtemp1) then
            fracsnow = 1.0
         else if (tair >= pxtemp2) then
             fracsnow = 0.0
         else

             ! Linear interpolate between 0 and 1
             fracsnow  = (tair - pxtemp1)/(pxtemp2 - pxtemp1)

             ! Make sure that fracsnow is between 0 and 1
             if (fracsnow > 1.) then
                 fracsnow = 1.
             else if (fracsnow < 0.) then
                 fracsnow = 0.
             end if

         end if
    end select

   ! now move on...
   fracrain = 1.0 - fracsnow

   ! snow phase of the precipitation
   pn = precip * fracsnow !* scf --  we are not using an undercatch factor (yet)

   ! rain phase of the precipitaiton
   rain = fracrain * precip

   ! wi is the accumulated swe stored in the ice portion of the snowpack
   ! (as opposed to liqiud water suspended)
   w_i =  w_i + pn

   ! Temperature and Heat deficit from new Snow
   ! The new snow temperature will be, at a minimum, zero
   if (tair < 0.0) then
       t_snow_new = tair
       ! delta_hd_snow = change in the heat deficit due to snowfall (mm)
       delta_hd_snow = - (t_snow_new * pn) / (80 / 0.5)
       t_rain = pxtemp
   else
       t_snow_new = 0.0
       delta_hd_snow = 0.0
       t_rain = tair
   end if


  ! Antecedent Temperature Index (ait)
  if (pn > (1.5 * dt)) then
      ait = t_snow_new
  else
      ait = ait + tipm_dt * (tair - ait)
  end if

  ! ensure that ait is always negative
  if (ait > 0) then
      ait = 0
  end if

  ! heat deficit
  ! this is the change in heat deficit caused by the temperature gradient (sensible heat...)
  ! *nmf* is an input parameter controlling this exchange
  ! notice that the melt factor (mf) also goes in here
  delta_hd_t = nmf * (dt / 6.0) * ((mf) / mfmax) * (ait - t_snow_new)


  ! Compute Rain-on-snow melt and Non-Rain Melt
  ! First, check if there is ice with which to melt...
  if (w_i > 0.) then
    !-----------------------------------------------------------------
    ! Rain-On-Snow Melt
    ! two options -- rain is intense enough for extra ROS melt or not
    !-----------------------------------------------------------------

    ! 1) if rain exceeds .25 mm/hour, then compute ROS melt
    ! notice that this function does not account for the amount
    ! of snow/heat content of the pack....
    if (rain  > (0.25 * dt)) then
      !melt (mm) during rain-on-snow periods is:
      m_ros1 = MAX(stefan * dt * (((tair + 273.0) ** 4) - (273.0 ** 4)), 0.0)
      m_ros2 = MAX(0.0125 * rain * t_rain, 0.0)
      m_ros3 = MAX(8.5 * uadj * (dt / 6.0) * ((0.9 * e_sat) - 6.11 + (0.00057 * p_atm * tair)), 0.0)
      m_ros = m_ros1 + m_ros2 + m_ros3

    ! 2) no rain on snow melt if the rain isn't intense enough
    else
      m_ros = 0.0
    end if

    !-----------------------------------------------------------------
    ! Non-Rain Melt
    ! two options -- rain is intense enough for extra ROS melt or not
    !-----------------------------------------------------------------
    if (rain <= (0.25 * dt).and.(tair > mbase)) then
        ! melt during non-rain periods is:
        m_nr = (mf * (tair - mbase)) + (0.0125 * rain * t_rain)
    else
        m_nr = 0.0
    end if

  ! Set melt to zero if there's no snow
  else
    m_ros = 0.0
    m_nr = 0.0
  end if


  ! Ripeness of the snow cover
  melt = m_ros + m_nr
  if (melt <= 0) then
      melt = 0.0
  end if

  if (melt < w_i) then
      w_i = w_i - melt
  else
      melt = w_i + w_q
      w_i = 0.0
  end if

  !qw = liquid water available melted/rained at the snow surface (mm)
  qw = melt + rain

  ! w_qx = liquid water capacity (mm)
  w_qx = plwhc * w_i

  ! deficit = heat deficit (mm)
  deficit = deficit + delta_hd_snow + delta_hd_t

  ! Determine limits of the heat deficit
  if (deficit < 0) then
      deficit = 0.0
  else if (deficit > (0.33 * w_i)) then
      deficit = 0.33 * w_i
  end if

  ! define the heat deficit + liquid water storage capacity
  ! !!!! why does it look like this exactly? !!!!
  k = deficit * (1 + plwhc) + w_qx

  ! If there is liquid stored as ice... ie a snowpack exists
  if (w_i > 0.0) then

      ! Case 1) Snow cover is completely "Ripe"
      if ((qw + w_q) > k) then
          ! Then the snow is RIPE
          ! # Excess liquid water (mm)
          ex_lw = qw + w_q - w_qx - (deficit * (1 + plwhc))
          ! # fills liquid water capacity
          w_q = w_qx
          ! # w_i increases because water refreezes as heat deficit is
          ! # decreased
          w_i = w_i + deficit
          deficit = 0.0

      ! Case 2) The snow is NOT yet ripe, but ice is being melted
      else if ((qw >= deficit).and.((qw + w_q) <= k)) then
          ex_lw = 0.0
          w_q = w_q + qw - deficit
          ! # w_i increases because water refreezes as heat deficit is
          ! # decreased
          w_i = w_i + deficit
          deficit = 0.0

      ! Case 3) Snow is not yet ripe
      else
          ex_lw = 0.0  ! there is no excess liquid water
          ! w_i increases because water refreezes as heat deficit is decreased
          w_i = w_i + qw
          deficit = deficit - qw
      end if

      ! Update the SWE, which is the sum of liquid + ice phase water
      swe = w_i + w_q

  ! Otherwise, there is no snowpack! melt
  else
      ! all of the ice is gone... so there is also no liquid
      ! water stored in the pack....
      w_q = 0
      ex_lw = qw
      swe = 0
  end if

  ! deficit...
  if (deficit == 0) then
      ait = 0
  end if

  ! End of model execution
  ! model_swe = swe  ! total swe (mm) at this time step
  outflow = ex_lw

end subroutine snow17



subroutine snow17driver(ntimes, jdayVec, precipVec, tairVec, & ! INPUTS
                        nlayers, dz,  dt, rvs,               & ! INPUTS            !layer_areas
                        OPG_METHOD, &
                        pm,  &                                 ! precip parameter
                        opg, &                                 ! precip parameter
                        bias, &                                ! precip parameter
                        uadj, &                                ! parameters
                        mbase, &                               ! parameters
                        mfmax, &                               ! parameters
                        mfmin, &                               ! parameters
                        tipm, &                                ! parameters
                        nmf, &                                 ! parameters
                        plwhc, &                               ! parameters
                        pxtemp, &                              ! parameters
                        pxtemp1, &                             ! parameters
                        pxtemp2, &                             ! parameters
                        alpha,   &
                        ecutoff,   &
                        t_lapse, &                             ! temperature lapse rate
                        t_bias, &                              ! temperature bias *C
                        outflowVec, &                          ! OUTPUT
                        sweVec, &                              ! OUTPUT
                        rainVec, ptotVec, tempVec)                                !  OUTPUT

                        implicit none

                        ! INPUT forcings
                        integer, intent(in) :: ntimes                    ! Number of timesteps..
                        real, intent(in), dimension(ntimes) :: jdayVec   ! julian day
                        real, intent(in), dimension(ntimes) :: precipVec ! precipitation
                        real, intent(in), dimension(ntimes) :: tairVec   ! air temperature

                        ! INPUT paramters
                        integer, intent(in) :: nlayers                   ! Number of model layers
                        real, intent(in), dimension(nlayers) :: dz       ! vector of *mean* elevtaion differences for each nlayers
                        !real, intent(in), dimension(nlayers) :: layer_areas   ! vector of elevation band areas

                        ! input parameters -- precipitation adjustment
                        integer, intent(in) :: OPG_METHOD
                        real, intent(in) :: pm
                        real, intent(in) :: opg
                        real, intent(in) :: bias

                        ! INPUT paramters -- Snow 17
                        real, intent(in) :: dt
                        integer, intent(in) :: rvs
                        real, intent(in) :: uadj
                        real, intent(in) :: mbase
                        real, intent(in) :: mfmax
                        real, intent(in) :: mfmin
                        real, intent(in) :: tipm
                        real, intent(in) :: nmf
                        real, intent(in) :: plwhc
                        real, intent(in) :: pxtemp
                        real, intent(in) :: pxtemp1
                        real, intent(in) :: pxtemp2
                        real, intent(in) :: alpha
                        real, intent(in) :: ecutoff
                        real, intent(in) :: t_lapse
                        real, intent(in) :: t_bias

                        ! OUTPUTS
                        real, intent(out), dimension(nlayers,ntimes) :: outflowVec
                        real, intent(out), dimension(nlayers,ntimes) :: sweVec
                        real, intent(out), dimension(nlayers,ntimes) :: rainVec
                        real, intent(out), dimension(nlayers,ntimes) :: ptotVec
                        real, intent(out), dimension(nlayers,ntimes) :: tempVec


                        !internal

                        ! states
                        !real, parameter :: t_lapse = -.0065    ! 6.5 c/m
                        real, parameter :: stat_elev = 30.9    ! elevation in HUNDREDS of meters. Approx elevation of the Butte Snotel
                        real, parameter :: opg_lower_adj = .25  ! reduced the opg gradient for lower stations
                        real  :: swe
                        real  :: ait
                        real  :: w_q
                        real  :: w_i
                        real  :: deficit
                        real  :: outflow
                        real ::  elevation
                        real ::  precip_il ! precipitation in the lth layer
                        real ::  temp_il
                        real ::  opg_adjusted
                        real ::  dz_adj
                        real ::  rain

                        ! misc
                        integer :: i   ! for time loop
                        integer :: l   ! for layer loop

                        ! loop through layers
                        do l=1,nlayers

                          ! Set Initial conditions for the layer
                          ait = 0.
                          swe = 0.
                          ait = 0.
                          w_q = 0.
                          w_i = 0.
                          deficit = 0.
                          outflow = 0.

                          ! midpoint elevation ...
                          elevation = stat_elev + dz(l)/100 ! elevation is in units of 100s of meters for whatever reason...

                          ! loop through timesteps
                          do i=1,ntimes

                            ! Adjust the precipitation
                            select case(OPG_METHOD)

                            ! basic linear adjustment
                            case(0)
                              ! adjust the precipitaiton
                              if (precipVec(i) > 0.) then
                                ! make sure that a negative bias doesn't cause negative precipitation
                                !precip_il = MAX(0., precipVec(i) + opg * dz(l) * precipVec(i) + bias)
                                precip_il = MAX(0., (10.**pm)*precipVec(i)*(1 +  opg * dz(l)) + bias)

                              else
                                precip_il = 0
                              end if

                            ! some kind of curvy adjustment
                            case(1)

                              ! adjust the precipitaiton IF there's precip
                              if (precipVec(i) > 0.) then

                                ! we are below the station. reduce the precip lapse rate ...
                                if (dz(l) < 0) then
                                  opg_adjusted = opg*opg_lower_adj
                                ! We are above the station. keep the precip lapse rate as it is
                                else
                                  opg_adjusted = opg
                                end if

                                !precip_il = MAX(0., precipVec(i) + opg_adjusted * dz(l) * precipVec(i) + bias)
                                precip_il = MAX(0., (10.**pm)*precipVec(i)*(1 +  opg_adjusted * dz(l)) + bias)

                              ! Keep precip zero if there's none
                              else
                                precip_il = 0
                              end if
                            !
                            case(2) !!! CUT OFF THE PRECIP ENHANCMENT UP HIGH -- above 3500 m (10s of meters in code)

                              ! adjust the precipitaiton
                              if (precipVec(i) > 0.) then

                                if (elevation > ecutoff)  then   ! this is 3500 m -- based on looking at ASO swe profile
                                  dz_adj = dz(l) - (elevation-ecutoff)*100*alpha
                                  precip_il = MAX(0., (10.**pm)*precipVec(i)*(1 +  opg * dz_adj) + bias)
                                 ! We are above the station. keep the precip lapse rate as it is
                                else
                                  ! make sure that a negative bias doesn't cause negative precipitation
                                  !precip_il = MAX(0., precipVec(i) + opg * dz(l) * precipVec(i) + bias)
                                  precip_il = MAX(0., (10.**pm)*precipVec(i)*(1 +  opg * dz(l)) + bias)
                                end if


                              else
                                precip_il = 0
                              end if

                            end select

                            ! Adjust air temperature
                            temp_il = tairVec(i) + t_lapse * dz(l) + t_bias

                            call snow17(jdayVec(i), &
                                        precip_il, &
                                        temp_il, &
                                        elevation, &
                                        dt, &
                                        rvs, &
                                        uadj, &
                                        mbase, &
                                        mfmax, &
                                        mfmin, &
                                        tipm, &
                                        nmf, &
                                        plwhc, &
                                        pxtemp, &
                                        pxtemp1, &
                                        pxtemp2, &
                                        swe, &
                                        ait, &
                                        w_q, &
                                        w_i, &
                                        deficit, &
                                        outflow, &  ! out
                                        rain)       ! out

                            ! store the data...
                            outflowVec(l,i) = outflow
                            sweVec(l,i) = swe
                            rainVec(l,i) = rain
                            ptotVec(l,i) = precip_il
                            tempVec(l,i) = temp_il
                        end do
                      end do


end subroutine snow17driver



end module snowmodule17





