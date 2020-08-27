! This module derived from the NWM repo found here: 
! https://github.com/NCAR/wrf_hydro_nwm_public/blob
! /master/trunk/NDHMS/Routing/Reservoirs/Level_Pool/module_levelpool.F

! This module defines and instantiates objects
! for a level pool type reservoir. The level
! pool reservoir type inherits input and
! output types from the reservoir base
! module and calls instantiation of these into
! sub-objects. The level pool reservoir type
! also points to types for level pool properties
! and state and calls instantiation of these into
! sub-objects. This module also contains the
! subroutine to run level pool reservoir that is
! derived from the reservoir base type interface
! to run reservoir. Running level pool will
! then call the LEVELPOOL_PHYSICS subroutine, which
! processes the given inputs, properties, and
! state for a particular level pool reservoir and
! returns the output/outflow.

module module_levelpool

    implicit none

contains

    ! ------------------------------------------------
    !   SUBROUTINE LEVELPOOL
    ! ------------------------------------------------

    subroutine levelpool_physics(dt,qi0,qi1,ql,ar,we,maxh,wc,wl,dl,oe,oc,oa,H0,H1,qo1)

        !! ----------------------------  argument variables
        !! All elevations should be relative to a common base (often belev(k))

        ! integer, intent(IN) :: ln      ! lake number
        real, intent(IN)    :: dt      ! routing period [s]
        real, intent(IN)    :: qi0     ! inflow at previous timestep (cms)
        real, intent(IN)    :: qi1     ! inflow at current timestep (cms)
        real, intent(IN)    :: ql      ! lateral inflow
        real, intent(IN)    :: ar      ! area of reservoir (km^2)
        real, intent(IN)    :: we      ! bottom of weir elevation
        real, intent(IN)    :: wc      ! weir coeff.
        real, intent(IN)    :: wl      ! weir length (m)
        real, intent(IN)    :: dl      ! dam length(m)
        real, intent(IN)    :: oe      ! orifice elevation
        real, intent(IN)    :: oc      ! orifice coeff.
        real, intent(IN)    :: oa      ! orifice area (m^2)
        real, intent(IN)    :: maxh    ! max depth of reservoir before overtop (m)
        real, intent(IN)    :: H0      ! water elevation height (m)
        real, intent(OUT)   :: H1      ! water elevation height (m)
        real, intent(OUT)   :: qo1     ! outflow at current timestep

        !!DJG Add lake option switch here...move up to namelist in future versions...
        integer :: LAKE_OPT            ! Lake model option (move to namelist later)
        real    :: H, Htmp             ! Temporary assign of incoming lake el. (m)

        !! ----------------------------  local variables
        real :: sap                    ! local surface area values
        real :: discharge              ! storage discharge m^3/s
        real :: tmp1, tmp2
        real :: dh, dh1, dh2, dh3      ! Depth in weir, and height function for 3 order RK
        real :: It, Itdt_3, Itdt_2_3   ! inflow hydrographs
        real :: maxWeirDepth           !maximum capacity of weir
        !real :: hdiff_vol, qdiff_vol   ! water balance check variables
        !! ----------------------------  subroutine body: from chow, mad mays. pg. 252
        !! -- determine from inflow hydrograph


        !!DJG Set hardwire for LAKE_OPT...move specification of this to namelist in
        !future versions...
        LAKE_OPT = 2
        Htmp = H0  !temporary set of incoming lake water elevation...
        !hdiff_vol = 0.0
        !qdiff_vol = 0.0

        !!DJG IF-block for lake model option  1 - outflow=inflow, 2 - Chow et al level
        !pool, .....
        if (LAKE_OPT == 1) then     ! If-block for simple pass through scheme....

           qo1 = qi1                 ! Set outflow equal to inflow at current time
           H = Htmp                  ! Set new lake water elevation to incoming lake el.

        else if (LAKE_OPT == 2) then   ! If-block for Chow et al level pool scheme

           It = qi0
           Itdt_3   = qi0 + ((qi1 + ql - qi0) * 0.33)
           Itdt_2_3 = qi0 + ((qi1 + ql - qi0) * 0.67)
           maxWeirDepth =  maxh - we

           !assume vertically walled reservoir
           !remove this when moving to a variable head area volume
           sap = ar * 1.0E6

           !-- determine Q(dh) from elevation-discharge relationship
           !-- and dh1
           dh = H - we
           if (dh > maxWeirDepth) then
              dh = maxWeirDepth
           endif

           tmp1 = oc * oa * sqrt(2. * 9.81 * ( H - oe )) !orifice at capacity
           tmp2 = wc * wl * (dh ** (3./2.))  !weir flows at capacity

           !determine the discharge based on current height
           if(H > maxh) then
             discharge =  tmp1 + tmp2 + (wc* (wl*dl) * (H-maxh)**(3./2.)) !overtop
           else if (dh > 0.0 ) then              !! orifice and weir discharge
             discharge = tmp1 + tmp2
           else if ( H > oe ) then     !! only orifice flow
             discharge = oc * oa * sqrt(2. * 9.81 * ( H - oe ) )
           else
             discharge = 0.0   !in the dead pool
           endif

           if (sap > 0) then
              dh1 = ((It - discharge)/sap)*dt
           else
              dh1 = 0.0
           endif

           !-- determine Q(H + dh1/3) from elevation-discharge relationship
           !-- dh2
           dh = (H+dh1/3) - we
           if (dh > maxWeirDepth) then
              dh = maxWeirDepth
           endif

           tmp1 = oc * oa * sqrt(2. * 9.81 * ( (H+dh1/3.) - oe ) )
           tmp2 = wc * wl * (dh ** (3./2.))

           !determine the discharge based on current height
           if(H > maxh) then
             discharge =  tmp1 + tmp2 + (wc* (wl*dl) * (H-maxh)**(3./2.)) !overtop
           else if (dh > 0.0 ) then              !! orifice and weir discharge
             discharge = tmp1 + tmp2
           else if ( (H+dh1/3) > oe ) then     !! only orifice flow,not full
             discharge = oc * oa * sqrt(2. * 9.81 * ( (H+dh1/3.) - oe ) )
           else
             discharge = 0.0
            endif


           if (sap > 0.0) then
              dh2 = ((Itdt_3 - discharge)/sap)*dt
           else
              dh2 = 0.0
           endif

           !-- determine Q(H + 2/3 dh2) from elevation-discharge relationship
           !-- dh3
           dh = (H + (0.667*dh2)) - we
           if (dh > maxWeirDepth) then
              dh = maxWeirDepth
           endif

           tmp1 = oc * oa * sqrt(2. * 9.81 * ( (H+dh2*0.667) - oe ) )
           tmp2 = wc * wl * (dh ** (3./2.))

           !determine the discharge based on current height
           if(H > maxh) then  ! overtop condition, not good!
              discharge =  tmp1 + tmp2 + (wc* (wl*dl) * (H-maxh)**(3./2.)) !overtop
           else if (dh > 0.0 ) then              !! orifice and weir discharge
              discharge = tmp1 + tmp2
           else if ( (H+dh2*0.667) > oe ) then     !! only orifice flow,not full
              discharge = oc * oa * sqrt(2. * 9.81 * ( (H+dh2*0.667) - oe ) )
           else
              discharge = 0.0
           endif

           if (sap > 0.0) then
              dh3 = ((Itdt_2_3 - discharge)/sap)*dt
           else
              dh3 = 0.0
           endif

           !-- determine dh and H
           dh = (dh1/4.) + (0.75*dh3)
           H = H + dh

           !-- compute final discharge
           dh = H - we
           if (dh > maxWeirDepth) then
              dh = maxWeirDepth
           endif

           tmp1 = oc * oa * sqrt(2. * 9.81 * ( H - oe ) )
           tmp2 = wc * wl * (dh ** (3./2.))

           !determine the discharge based on current height
           if(H > maxh) then  ! overtop condition, not good!
              discharge =  tmp1 + tmp2 + (wc* (wl*dl) * (H-maxh)**(3./2.)) !overtop
           else if (dh > 0.0 ) then              !! orifice and overtop discharge
              discharge = tmp1 + tmp2
           else if ( H > oe ) then     !! only orifice flow,not full
              discharge = oc * oa * sqrt(2. * 9.81 * ( H - oe ) )
           else
              discharge = 0.0
           endif

           H1 = H
           qo1  = discharge  ! return the flow rate from reservoir

        !#ifdef HYDRO_D
        !#ifndef NCEP_WCOSS
        !   ! Water balance check
        !   qdiff_vol = (qi1+ql-qo1)*dt !m3
        !   hdiff_vol = (H-Htmp)*sap    !m3
        !22 format(f8.4,2x,f8.4,2x,f8.4,2x,f8.4,2x,f8.4,2x,f6.0,2x,f20.1,2x,f20.1)
        !   open (unit=67, &
        !     file='lake_massbalance_out.txt', status='unknown',position='append')
        !   write(67,22) Htmp, H, qi1, ql, qo1, dt, qdiff_vol, hdiff_vol
        !   close(67)
        !#endif
        !#endif

        !23 format('botof H dh orf wr Q',f8.4,2x,f8.4,2x,f8.3,2x,f8.3,2x,f8.2)
        !24 format('ofonl H dh sap Q ',f8.4,2x,f8.4,2x,f8.0,2x,f8.2)


        else   ! ELSE for LAKE_OPT....
        endif  ! ENDIF for LAKE_OPT....

        return

    ! ----------------------------------------------------------------
    end subroutine levelpool_physics
    ! ----------------------------------------------------------------

end module module_levelpool
