module submuskingcunge_wrf_module

    use precis
    implicit none

contains

subroutine submuskingcunge(    &
     qup,  quc, &
     qdp, ql,   dt,  So,   dx, &
     n,   Cs,   Bw,  Tw, TwCC, &
     nCC, depthp, qdc, velc, depthc)

        IMPLICIT NONE

        REAL, intent(IN)       :: dt         ! routing period in  seconds
        REAL, intent(IN)       :: qup        ! flow upstream previous timestep
        REAL, intent(IN)       :: quc        ! flow upstream current timestep
        REAL, intent(IN)       :: qdp        ! flow downstream previous timestep
        REAL, intent(OUT)      :: qdc        ! flow downstream current timestep
        REAL, intent(IN)       :: ql         ! lateral inflow through reach (m^3/sec)
        REAL, intent(IN)       :: Bw         ! bottom width (meters)
        REAL, intent(IN)       :: Tw         ! top width before bankfull (meters)
        REAL, intent(IN)       :: TwCC       ! top width of Compund (meters)
        REAL, intent(IN)       :: nCC        ! mannings of compund
        REAL, intent(IN)       :: Cs         ! Channel side slope slope
        REAL, intent(IN)       :: So         ! Channel bottom slope %
        REAL, intent(IN)       :: dx         ! channel lngth (m)
        REAL, intent(IN)       :: n          ! mannings coefficient
        REAL, intent(OUT)      :: velc       ! channel velocity
        REAL, intent(IN)       :: depthp     ! depth of flow in channel
        REAL, intent(OUT)      :: depthc     ! depth of flow in channel

!--local variables
        REAL    :: C1, C2, C3, C4
        REAL    :: Km             !K travel time in hrs in reach
        REAL    :: X              !weighting factors 0<=X<=0.5
        REAL    :: Ck             ! wave celerity (m/s)

!-- channel geometry and characteristics, local variables
        REAL    :: Twl            ! top width at simulated flow (m)
        REAL    :: AREA,AREAC     ! Cross sectional area channel and compound m^2
        REAL    :: Z              ! trapezoid distance (m)
        REAL    :: R,RC           ! Hydraulic radius of channel and compound
        REAL    :: WP,WPC         ! wetted perimmeter of channel and compound
        REAL    :: h              ! depth of flow in channel
        REAL    :: h_0,h_1        ! secant method estimates
        REAL    :: bfd            ! bankfull depth (m)
        REAL    :: Qj_0           ! secant method estimates
        REAL    :: Qj             ! intermediate flow estimate
        REAL    :: D,D1           ! diffusion coeff
        REAL    :: dtr            ! required timestep, minutes
        REAL    :: aerror,rerror  ! absolute and relative error
        REAL    :: hp             ! courant, previous height
        INTEGER :: iter, maxiter  ! maximum number of iterations

!-- local variables.. needed if channel is sub-divded
        REAL    :: a,b,c,F
        REAL    :: mindepth   !  minimum depth in channel
        INTEGER :: i,tries    !-- channel segment counter

        maxiter  = 100
        mindepth = 0.01

        aerror = 0.01
        rerror = 1.0
        tries = 0

!-------------  locals
        if(Cs .eq. 0.00000000) then
         z = 1.0
        else
         z = 1.0/Cs          !channel side distance (m)
        endif

        if(Bw .gt. Tw) then   !effectively infinite deep bankful
           bfd = Bw/0.00001
        elseif (Bw .eq. Tw) then
          bfd =  Bw/(2.0*z)  !bankfull depth is effectively
        else
          bfd =  (Tw - Bw)/(2.0*z)  !bankfull depth (m)
        endif
        !print *,"bfd:",bfd
        !qC = quc + ql !current upstream in reach

        if (n .le. 0.0 .or. So .le. 0.0 .or. z .le. 0.0 .or. Bw .le. 0.0) then
          print*, "Error in channel coefficients -> Muskingum cunge", n, So, z, Bw
          !call hydro_stop("In MUSKINGCUNGE() - Error in channel coefficients")
        end if

!-------------  Secant Method
        depthc = max(depthp, 0.0)
        h     = (depthc * 1.33) + mindepth !1.50 of  depth
        h_0   = (depthc * 0.67)            !0.50 of depth

        if(ql .gt. 0.0 .or. qup .gt. 0.0 .or. qdp .gt. 0.0) then  !only solve if there's water to flux

110     continue

        Qj_0  = 0.0                       !- initial flow of lower interval
        !WPC    = 0.0_prec
        !AREAC  = 0.0_prec
        iter   = 0

        do while (rerror .gt. 0.01 .and. aerror .ge. mindepth .and. iter .le. maxiter)

           AREAC  = 0.0_prec
           WPC    = 0.0_prec

          !----- lower interval  --------------------
           Twl = Bw + 2.0*z*h_0      !--top surface water width of the channel inflow

            if(h_0 .gt. bfd) then !water outside of defined channel
             AREA =  (Bw + bfd * z) * bfd
             AREAC = (TwCC * (h_0 -bfd)) !assume compound component is rect. chan, that's 3 times the Tw
             WP = (Bw + 2.0 * bfd * sqrt(1.0 + z*z))
             WPC = TwCC + (2.0 * (h_0-bfd)) !WPC is 2 times the Tw
             R   = (AREA + AREAC)/(WP +WPC)  ! hydraulic radius
            else
              AREA = (Bw + h_0 * z ) * h_0
              WP = (Bw + 2.0 * h_0 * sqrt(1.0 + z*z))

              if(WP .gt. 0.0) then
               R = AREA/ WP
              else
               R = 0.0
              endif

            endif

          if(h_0 .gt. bfd) then !water outside of defined channel
            !weight the celerity by the contributing area, and assume that the mannings
            !of the spills is 2x the manning of the channel
                Ck = max(0.0,((sqrt(So)/n)*((5./3.)*R**(2./3.) - &
                ((2./3.)*R**(5./3.)*(2.0*sqrt(1.0 + z*z)/(Bw+2.0*bfd*z))))*AREA &
                + ((sqrt(So)/(nCC))*(5./3.)*(h_0-bfd)**(2./3.))*AREAC)/(AREA+AREAC))
          else
               if(h_0 .gt. 0.0) then
                 Ck = max(0.0,(sqrt(So)/n)*((5./3.)*R**(2./3.) - &
                 ((2./3.)*R**(5./3.)*(2.0*sqrt(1.0 + z*z)/(Bw+2.0*h_0*z)))))
                else
                 Ck = 0.0
                endif
          endif

          if(Ck .gt. 0.000000) then
            Km = max(dt,dx/Ck)
          else
            Km = dt
          endif

          if(h_0 .gt. bfd) then !water outside of defined channel
             X = min(0.5,max(0.0,0.5*(1-(Qj_0/(2.0*TwCC*So*Ck*dx)))))
          else
            if(Ck .gt. 0.0) then
              X = min(0.5,max(0.0,0.5*(1-(Qj_0/(2.0*Twl*So*Ck*dx)))))
            else
              X = 0.5
            endif
          endif

           D = (Km*(1.000 - X) + dt/2.0000)              !--seconds
            if(D .eq. 0.0) then
              print *, "FATAL ERROR: D is 0 in MUSKINGCUNGE", Km, X, dt,D
              !call hydro_stop("In MUSKINGCUNGE() - D is 0.")
           endif

           C1 =  (Km*X + dt/2.000000)/D
           C2 =  (dt/2.0000 - Km*X)/D
           C3 =  (Km*(1.00000000-X)-dt/2.000000)/D
!          C1 =  max(0.0,min(1.0,1.0-C3))
           C4 =  (ql*dt)/D
!          C4 =  (ql*dt)/D - (ChannK * (WP + WPC) * dx)  !-- DY & LKR lat inflow + channel loss

           !!Uncomment to show WP/WPC behavior above bankfull
           !print *,"secant1 --", "WP:", WP, "WPC:", WPC
           if((WP+WPC) .gt. 0.0) then  !avoid divide by zero
             Qj_0 =  ((C1*qup)+(C2*quc)+(C3*qdp) + C4) - ((1/(((WP*n)+(WPC*nCC))/(WP+WPC))) * &
                    (AREA+AREAC) * (R**(2./3.)) * sqrt(So)) !f0(x)
           endif

           AREAC  = 0.0_prec
           WPC    = 0.0_prec

           !--upper interval -----------
           Twl = Bw + 2.0*z*h                    !--top width of the channel inflow

           if(h .gt. bfd) then !water outside of defined channel
             AREA =  (Bw + bfd * z) * bfd
             AREAC = (TwCC * (h-bfd)) !assume compound component is rect. chan, that's 3 times the Tw
             WP = (Bw + 2.0 * bfd * sqrt(1.0 + z*z))
             WPC = TwCC + (2.0*(h-bfd)) !the additional wetted perimeter
             R   = (AREA + AREAC)/(WP +WPC)
            ! RC  = AREAC/WPC
           else
              AREA = (Bw + h * z ) * h
              WP = (Bw + 2.0 * h * sqrt(1.000000 + z*z))
              if(WP .gt. 0.0) then
               R = AREA/WP
              else
               R = 0.0
              endif
           endif

          if(h .gt. bfd) then !water outside of defined channel, assumed rectangular, 3x TW and n = 3x
                Ck = max(0.0,((sqrt(So)/n)*((5./3.)*R**(2./3.) - &
                ((2./3.)*R**(5./3.)*(2.0*sqrt(1.0 + z*z)/(Bw + 2.0*bfd*z))))*AREA &
                + ((sqrt(So)/(nCC))*(5./3.)*(h-bfd)**(2./3.))*AREAC)/(AREA+AREAC))
          else
               if(h .gt. 0.0) then !avoid divide by zero
                 Ck = max(0.0,(sqrt(So)/n)*((5./3.)*R**(2./3.) - &
                 ((2./3.)*R**(5./3.)*(2.0 * sqrt(1.0 + z*z)/(Bw + 2.0*h*z)))))
               else
                 Ck = 0.0
               endif
          endif
           if(Ck .gt. 0.0) then
            Km = max(dt,dx/Ck)
           else
            Km = dt
           endif
          if(h .gt. bfd) then !water outside of defined channel
            X = min(0.5,max(0.25,0.5*(1-(((C1*qup)+(C2*quc)+(C3*qdp) + C4)/(2.0*TwCC*So*Ck*dx)))))

          else

            if(Ck .gt. 0.0) then
             X = min(0.5,max(0.25,0.5*(1-(((C1*qup)+(C2*quc)+(C3*qdp) + C4)/(2.0*Twl*So*Ck*dx)))))
            else
             X = 0.5
            endif

          endif

           D = (Km*(1 - X) + dt/2)              !--seconds
            if(D .eq. 0.0) then
              print *, "FATAL ERROR: D is 0 in MUSKINGCUNGE", Km, X, dt,D
              !call hydro_stop("In MUSKINGCUNGE() - D is 0.")
           endif

           C1 =  (Km*X + dt/2.000000)/D
           C2 =  (dt/2.000000 - Km*X)/D
           C3 =  (Km*(1.000000-X)-dt/2.000000)/D
!          C1 =  max(0.0,min(1.0,1.0-C3)) !! eliminate influence of upstream current
           C4 =  (ql*dt)/D
!          C4 =  (ql*dt)/D  - (ChannK * dx * (WP+WPC))  !-- (loss units: m/s * m * m)

           if( (C4 .lt. 0.0) .and. (abs(C4) .gt. (C1*qup)+(C2*quc)+(C3*qdp)))  then
            C4 = -((C1*qup)+(C2*quc)+(C3*qdp))
           endif

           !!Uncomment to show WP/WPC behavior above bankfull
           !print *,"secant2 --", "WP:", WP, "WPC:", WPC
           if((WP+WPC) .gt. 0.0) then
            Qj =  ((C1*qup)+(C2*quc)+(C3*qdp) + C4) - ((1.0000000/(((WP*n)+(WPC*nCC))/(WP+WPC))) * &
                    (AREA+AREAC) * (R**(2./3.)) * sqrt(So)) !f(x)
           endif

           if(Qj_0-Qj .ne. 0.0) then
             h_1 = h - ((Qj * (h_0 - h))/(Qj_0 - Qj)) !update h, 3rd estimate
              if(h_1 .lt. 0.0) then
                h_1 = h
              endif
           else
             h_1 = h
           endif

           if(h .gt. 0.0) then
             rerror = abs((h_1 - h)/h) !relative error is new estatimate and 2nd estimate
             aerror = abs(h_1 -h)      !absolute error
           else
             rerror = 0.0
             aerror = 0.9
           endif

           h_0  = max(0.0,h)
           h    = max(0.0,h_1)
           iter = iter + 1

           if( h .lt. mindepth) then  ! exit loop if depth is very small
             goto 111
           endif

         end do

111      continue

         if(iter .ge. maxiter) then

           tries = tries + 1
           if(tries .le. 4) then  ! expand the search space
             h     =  h * 1.33
             h_0   =  h_0 * 0.67
             maxiter = maxiter + 25 !and increase the number of allowable iterations
            goto 110
           endif

           print*, "Musk Cunge WARNING: Failure to converge"
           !print*, 'RouteLink index:', linkls_s(my_id+1) - 1
           print*, "id,err,iters,tries", rerror, iter, tries
           print*, "Ck,X,dt,Km",Ck,X,dt,Km
           print*, "So,dx,h",So,dx,h
           print*, "qup,quc,qdp,ql", qup,quc,qdp,ql
           print*, "bfd,Bw,Tw,Twl", bfd,Bw,Tw,Twl
           print*, "Qmc,Qmn", (C1*qup)+(C2*quc)+(C3*qdp) + C4,((1/(((WP*n)+(WPC*nCC))/(WP+WPC))) * &
                    (AREA+AREAC) * (R**(2./3.)) * sqrt(So))
         endif

!yw added for test
      !DY and LKR Added to update for channel loss
        if(((C1*qup)+(C2*quc)+(C3*qdp) + C4) .lt. 0.0) then
!       MUSKINGCUNGE =  MAX( ( (C1*qup)+(C2*quc) + C4),((C1*qup)+(C3*qdp) + C4) )
           if( (C4 .lt. 0.0) .and. (abs(C4) .gt. (C1*qup)+(C2*quc)+(C3*qdp)) )  then ! channel loss greater than water in chan
             qdc = 0.0
           else
             qdc = MAX( ( (C1*qup)+(C2*quc) + C4),((C1*qup)+(C3*qdp) + C4) )
           endif
        else
!       MUSKINGCUNGE =  ((C1*qup)+(C2*quc)+(C3*qdp) + C4) !-- pg 295 Bedient huber
          qdc =  ((C1*qup)+(C2*quc)+(C3*qdp) + C4) !-- pg 295 Bedient huber

        endif

        Twl = Bw + (2.0*z*h)
        R = (h*(Bw + Twl) / 2.0) / (Bw + 2.0*(((Twl - Bw) / 2.0)**2.0 + h**2)**0.5)
        velc =  (1./n) * (R **(2.0/3.0)) * sqrt(So)  ! average velocity in m/s
        depthc = h

      else   ! no flow to route
       qdc = 0.0
       depthc = 0.0
     endif

! ----------------------------------------------------------------
END SUBROUTINE SUBMUSKINGCUNGE
! ----------------------------------------------------------------

end module submuskingcunge_wrf_module
