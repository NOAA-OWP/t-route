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

    use precis
    implicit none

contains

    ! ------------------------------------------------
    !   SUBROUTINE LEVELPOOL
    ! ------------------------------------------------

    subroutine levelpool_physics(dt,qi0,qi1,ql,ar,we,maxh,wc,wl,dl,oe,oc,oa,H0,H1,qo1)

        !! ----------------------------  argument variables
        !! All elevations should be relative to a common base (often belev(k))

        ! integer, intent(IN) :: ln      ! lake number
        real(prec), intent(IN)    :: dt      ! routing period [s]
        real(prec), intent(IN)    :: qi0     ! inflow at previous timestep (cms)
        real(prec), intent(IN)    :: qi1     ! inflow at current timestep (cms)
        real(prec), intent(IN)    :: ql      ! lateral inflow
        real(prec), intent(IN)    :: ar      ! area of reservoir (km^2)
        real(prec), intent(IN)    :: we      ! bottom of weir elevation
        real(prec), intent(IN)    :: wc      ! weir coeff.
        real(prec), intent(IN)    :: wl      ! weir length (m)
        real(prec), intent(IN)    :: dl      ! dam length(m)
        real(prec), intent(IN)    :: oe      ! orifice elevation
        real(prec), intent(IN)    :: oc      ! orifice coeff.
        real(prec), intent(IN)    :: oa      ! orifice area (m^2)
        real(prec), intent(IN)    :: maxh    ! max depth of reservoir before overtop (m)
        real(prec), intent(IN)    :: H0      ! water elevation height (m)
        real(prec), intent(OUT)   :: H1      ! water elevation height (m)
        real(prec), intent(OUT)   :: qo1     ! outflow at current timestep

        !!DJG Add lake option switch here...move up to namelist in future versions...
        integer :: LAKE_OPT            ! Lake model option (move to namelist later)
        real(prec)    :: H, Htmp             ! Temporary assign of incoming lake el. (m)

        !! ----------------------------  local variables
        real(prec) :: sap                    ! local surface area values
        real(prec) :: discharge              ! storage discharge m^3/s
        real(prec) :: tmp1, tmp2
        real(prec) :: dh, dh1, dh2, dh3      ! Depth in weir, and height function for 3 order RK
        real(prec) :: It, Itdt_3, Itdt_2_3   ! inflow hydrographs
        real(prec) :: maxWeirDepth           !maximum capacity of weir
        real(prec) :: exp32, quot23, zero, gravity  ! internal variables
        !real :: hdiff_vol, qdiff_vol   ! water balance check variables
        !! ----------------------------  subroutine body: from chow, mad mays. pg. 252
        !! -- determine from inflow hydrograph


        !!DJG Set hardwire for LAKE_OPT...move specification of this to namelist in
        !future versions...
        LAKE_OPT = 2
        Htmp = H0  !temporary set of incoming lake water elevation...
        exp32 = 3._prec/2._prec
        quot23 = 0.667_prec
        zero = 0.0_prec
        gravity = 9.81_prec
        !hdiff_vol = zero
        !qdiff_vol = zero

        !!DJG IF-block for lake model option  1 - outflow=inflow, 2 - Chow et al level
        !pool, .....
        if (LAKE_OPT == 1) then     ! If-block for simple pass through scheme....

           qo1 = qi1                 ! Set outflow equal to inflow at current time
           H = Htmp                  ! Set new lake water elevation to incoming lake el.

        else if (LAKE_OPT == 2) then   ! If-block for Chow et al level pool scheme

           It = qi0
           Itdt_3   = qi0 + ((qi1 + ql - qi0) * 0.33_prec)
           Itdt_2_3 = qi0 + ((qi1 + ql - qi0) * 0.67_prec)
           maxWeirDepth =  maxh - we

           !assume vertically walled reservoir
           !remove this when moving to a variable head area volume
           sap = ar * 1.0E6_prec

           !-- determine Q(dh) from elevation-discharge relationship
           !-- and dh1
           dh = H - we
           if (dh > maxWeirDepth) then
              dh = maxWeirDepth
           endif

           tmp1 = oc * oa * sqrt(2._prec * gravity * ( H - oe )) !orifice at capacity
           tmp2 = wc * wl * (dh ** (exp32))  !weir flows at capacity

           !determine the discharge based on current height
           if(H > maxh) then
             discharge =  tmp1 + tmp2 + (wc* (wl*dl) * (H-maxh)**(exp32)) !overtop
           else if (dh > zero) then              !! orifice and weir discharge
             discharge = tmp1 + tmp2
           else if ( H > oe ) then     !! only orifice flow
             discharge = oc * oa * sqrt(2._prec * gravity * ( H - oe ) )
           else
             discharge = zero   !in the dead pool
           endif

           if (sap > zero) then
              dh1 = ((It - discharge)/sap)*dt
           else
              dh1 = zero
           endif

           !-- determine Q(H + dh1/3) from elevation-discharge relationship
           !-- dh2
           dh = (H+dh1/3._prec) - we
           if (dh > maxWeirDepth) then
              dh = maxWeirDepth
           endif

           tmp1 = oc * oa * sqrt(2._prec * gravity * ( (H+dh1/3._prec) - oe ) )
           tmp2 = wc * wl * (dh ** (exp32))

           !determine the discharge based on current height
           if(H > maxh) then
             discharge =  tmp1 + tmp2 + (wc* (wl*dl) * (H-maxh)**(exp32)) !overtop
           else if (dh > zero) then              !! orifice and weir discharge
             discharge = tmp1 + tmp2
           else if ( (H+dh1/3._prec) > oe ) then     !! only orifice flow,not full
             discharge = oc * oa * sqrt(2._prec * gravity * ( (H+dh1/3._prec) - oe ) )
           else
             discharge = zero
            endif


           if (sap > zero) then
              dh2 = ((Itdt_3 - discharge)/sap)*dt
           else
              dh2 = zero
           endif

           !-- determine Q(H + 2/3 dh2) from elevation-discharge relationship
           !-- dh3
           dh = (H + (quot23*dh2)) - we
           if (dh > maxWeirDepth) then
              dh = maxWeirDepth
           endif

           tmp1 = oc * oa * sqrt(2._prec * gravity * ( (H+dh2*quot23) - oe ) )
           tmp2 = wc * wl * (dh ** (exp32))

           !determine the discharge based on current height
           if(H > maxh) then  ! overtop condition, not good!
              discharge =  tmp1 + tmp2 + (wc* (wl*dl) * (H-maxh)**(exp32)) !overtop
           else if (dh > zero) then              !! orifice and weir discharge
              discharge = tmp1 + tmp2
           else if ( (H+dh2*quot23) > oe ) then     !! only orifice flow,not full
              discharge = oc * oa * sqrt(2._prec * gravity * ( (H+dh2*quot23) - oe ) )
           else
              discharge = zero
           endif

           if (sap > zero) then
              dh3 = ((Itdt_2_3 - discharge)/sap)*dt
           else
              dh3 = zero
           endif

           !-- determine dh and H
           dh = (dh1/4._prec) + (0.75_prec*dh3)
           H = H + dh

           !-- compute final discharge
           dh = H - we
           if (dh > maxWeirDepth) then
              dh = maxWeirDepth
           endif

           tmp1 = oc * oa * sqrt(2._prec * gravity * ( H - oe ) )
           tmp2 = wc * wl * (dh ** (exp32))

           !determine the discharge based on current height
           if(H > maxh) then  ! overtop condition, not good!
              discharge =  tmp1 + tmp2 + (wc* (wl*dl) * (H-maxh)**(exp32)) !overtop
           else if (dh > zero) then              !! orifice and overtop discharge
              discharge = tmp1 + tmp2
           else if ( H > oe ) then     !! only orifice flow,not full
              discharge = oc * oa * sqrt(2._prec * gravity * ( H - oe ) )
           else
              discharge = zero
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
