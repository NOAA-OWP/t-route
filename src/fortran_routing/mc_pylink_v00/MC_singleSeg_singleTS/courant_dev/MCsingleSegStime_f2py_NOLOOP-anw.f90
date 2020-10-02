module muskingcunge_module

    ! the 'precis' module is build in varPrecision.f90
    ! sets the precision of real variables
    use precis
    implicit none

contains

subroutine muskingcungenwm(dt, qup, quc, qdp, ql, dx, bw, tw, twcc,&
    n, ncc, cs, s0, velp, depthp, qdc, velc, depthc, ck, cn, X)

    !* exactly follows SUBMUSKINGCUNGE in NWM:
    !* 1) qup and quc for a reach in upstream limit take zero values all the time
    !* 2) initial value of depth of time t of each reach is equal to the value at time t-1
    !* 3) qup as well as quc at time t for a downstream reach in a serial network takes
    !*    exactly the same value qdp at time t (or qdc at time t-1) for the upstream reach

    implicit none

    real(prec), intent(in) :: dt 
    real(prec), intent(in) :: qup, quc, qdp, ql
    real(prec), intent(in) :: dx, bw, tw, twcc, n, ncc, cs, s0
    real(prec), intent(in) :: velp
    real(prec), intent(in) :: depthp
    real(prec), intent(out) :: qdc, velc, depthc
    real(prec), intent(out) :: ck, cn, X
    real(prec) :: z
    real(prec) :: bfd, C1, C2, C3, C4

    !Uncomment next line for old initialization
    !real(prec) :: WPC, AREAC

    integer :: iter
    integer :: maxiter, tries
    real(prec) :: mindepth, aerror, rerror
    real(prec) :: R, twl, h_1, h, h_0, Qj, Qj_0

    ! qdc = 0.0
    ! velc = velp
    ! depthc = depthp

    !* parameters of Secant method
    maxiter  = 100
    mindepth = 0.01_prec

    aerror = 0.01_prec
    rerror = 1.0_prec
    tries = 0

    ! calculate side distance (m)
    if(cs .eq. 0.0_prec) then
        z = 1.0_prec
    else
        z = 1.0_prec/cs          !channel side distance (m)
    endif

    ! calculate bankfull depth
    if(bw .gt. tw) then   !effectively infinite deep bankful
        bfd = bw/0.00001_prec
    elseif (bw .eq. tw) then
        bfd =  bw/(2.0_prec*z)  !bankfull depth is effectively
    else
        bfd =  (tw - bw)/(2.0_prec*z)  !bankfull depth (m)
    endif

    ! Throw error if parameters are less than zero
    if (n .le. 0.0_prec .or. s0 .le. 0.0_prec .or. z .le. 0.0_prec .or. bw .le. 0.0_prec) then
        !print*, "Error in channel coefficients -> Muskingum cunge", n, s0, z, bw
        !call hydro_stop("In MUSKINGCUNGE() - Error in channel coefficients")
    end if

    ! ?? what is the diff between depthc and depthp ??
    depthc = max(depthp, 0.0_prec)
    
    ! initialize the value of h and h_0
    h     = (depthc * 1.33_prec) + mindepth ! depth of flow in channel
    h_0   = (depthc * 0.67_prec)            ! secant method estimate

    ! only solve if there's water to flux
    if(ql .gt. 0.0_prec .or. qup .gt. 0.0_prec .or. qdp .gt. 0.0_prec .or. qdc .gt. 0.0_prec) then  
110 continue

        !Uncomment next two lines for old initialization
        !WPC = 0.0_prec
        !AREAC = 0.0_prec

        iter = 0

        ! itteratively solve for depth ??
        do while (rerror .gt. 0.01_prec .and. aerror .ge. mindepth .and. iter .le. maxiter)

            !Uncomment next four lines for old initialization
            !call secant2_h(z, bw, bfd, twcc, s0, n, ncc, dt, dx, &
            !    qdp, ql, qup, quc, h_0, 1, WPC, Qj_0, C1, C2, C3, C4)
            !call secant2_h(z, bw, bfd, twcc, s0, n, ncc, dt, dx, &
            !    qdp, ql, qup, quc, h, 2, WPC, Qj, C1, C2, C3, C4)

            ! interval = 1 & h = h_0, return Qj_0, C1, C2, C3, C4
            call secant2_h(z, bw, bfd, twcc, s0, n, ncc, dt, dx, &
                qdp, ql, qup, quc, h_0, 1, Qj_0, C1, C2, C3, C4, X)
                
            ! interval = 2 & h = h, return Qj, C1, C2, C3, C4  
            call secant2_h(z, bw, bfd, twcc, s0, n, ncc, dt, dx, &
                qdp, ql, qup, quc, h, 2, Qj, C1, C2, C3, C4, X)

            if(Qj_0-Qj .ne. 0.0_prec) then
                h_1 = h - ((Qj * (h_0 - h))/(Qj_0 - Qj)) !update h, 3rd estimate

                if(h_1 .lt. 0.0_prec) then
                    h_1 = h
                endif
            else
                h_1 = h
            endif

            if(h .gt. 0.0_prec) then
                rerror = abs((h_1 - h)/h) !relative error is new estimate and 2nd estimate
                aerror = abs(h_1 -h)      !absolute error
            else
                rerror = 0.0_prec
                aerror = 0.9_prec
            endif

            ! update depths for secant method
            h_0  = max(0.0_prec,h)
            h    = max(0.0_prec,h_1)
            iter = iter + 1

            if( h .lt. mindepth) then  ! exit loop if depth is very small
                goto 111
            endif
            
        end do !*do while (rerror .gt. 0.01 .and. ....
        
111    continue

        ! if maximum iterations are exceeded, then try again
        if(iter .ge. maxiter) then
            tries = tries + 1

            ! expand the search space and increase the number of allowable itterations
            if(tries .le. 4) then  
                h     =  h * 1.33_prec
                h_0   =  h_0 * 0.67_prec
                maxiter = maxiter + 25
                goto 110
                
            endif
            
        endif


        ! calculate channel outflow, current timestep (i.e. t+1)
        if(((C1*qup)+(C2*quc)+(C3*qdp) + C4) .lt. 0.0_prec) then
            if( (C4 .lt. 0.0_prec) .and. (abs(C4) .gt. (C1*qup)+(C2*quc)+(C3*qdp)) )  then ! channel loss greater than water in chan
                qdc = 0.0_prec

            else
                qdc = MAX( ( (C1*qup)+(C2*quc) + C4),((C1*qup)+(C3*qdp) + C4) )

            endif
        else
            ! outflow(t+1) = [C1*inflow(t)] + [C2*inflow(t+1)] + [C3*outflow(t)] + C4
            qdc = ((C1*qup)+(C2*quc)+(C3*qdp) + C4) !-- pg 295 Bedient huber

        endif

        ! Hydraulic radius calculation for trapezoidal channel
        ! ??? Why is the velocity calculation based on a trapezoidal channel shape, while the depth and flow calcs assume compound channel?
        twl = bw + (2.0_prec*z*h)
        R = (h*(bw + twl) / 2.0_prec) / (bw + 2.0_prec*(((twl - bw) / 2.0_prec)**2.0_prec + h**2.0_prec)**0.5_prec)
        
        ! velocity calculation, via Manning's eqn.
        velc = (1.0_prec/n) * (R **(2.0_prec/3.0_prec)) * sqrt(s0)  !*average velocity in m/s
        
        ! depth at  the current timestep (t+1) is equal to the itteratively solved for depth, h
        depthc = h
        
    else   ! no flow to route
        qdc = 0.0_prec

        depthc = 0.0_prec
        
    end if !*if(ql .gt. 0.0 .or. ...
    
    ! *************************************************************
    ! call courant subroutine here
    ! *************************************************************
    call courant(h, bfd, bw, twcc, ncc, s0, n, z, dx, dt, ck, cn)

end subroutine muskingcungenwm

!**---------------------------------------------------**!
!*                                                     *!
!*                 SECANT2 SUBROUTINE                  *!
!*                                                     *!
!**---------------------------------------------------**!
!Uncomment this function signature for old initialization
!subroutine secant2_h(z, bw, bfd, twcc, s0, n, ncc, dt, dx, &
!    qdp, ql, qup, quc, h, interval, WPC, Qj, C1, C2, C3, C4)

!Uncomment this function signature for new initialization 
subroutine secant2_h(z, bw, bfd, twcc, s0, n, ncc, dt, dx, &
    qdp, ql, qup, quc, h, interval, Qj, C1, C2, C3, C4, X)

    implicit none

    real(prec), intent(in) :: z, bw, bfd, twcc, s0, n, ncc
    real(prec), intent(in) :: dt, dx
    real(prec), intent(in) :: qdp, ql, qup, quc
    real(prec), intent(in) :: h
    real(prec), intent(out) :: Qj, C1, C2, C3, C4, X
    integer,    intent(in) :: interval

    real(prec) :: twl, AREA, WP, R, Ck, Km, D
    integer    :: upper_interval, lower_interval

    !Uncomment for old initialization
    !real(prec), intent(out) :: WPC
    !real(prec) :: AREAC
    !Uncomment for new initialization 
    real(prec) :: WPC, AREAC

    twl = 0.0_prec
    WP = 0.0_prec

    !Uncomment next line for old initialization
    !AREA = 0.0_prec
    !Uncomment next two lines for new initialization 
    WPC = 0.0_prec
    AREAC = 0.0_prec

    R = 0.0_prec
    Ck = 0.0_prec
    Km = 0.0_prec
    X = 0.0_prec
    D = 0.0_prec

    !--upper interval -----------
    upper_interval = 1
    !--lower interval -----------
    lower_interval = 2

    ! top surface water width of the channel inflow
    twl = bw + 2.0_prec*z*h

    !**hydraulic radius, R
    
    ! if depth is greater than bankful, then water is outside of defined channel
    if(h .gt. bfd) then
    
        ! channel area, main channel, trapezoid
        AREA =  (bw + bfd * z) * bfd
        
        ! channel area, compound (floodplain), rectangular
        AREAC = (twcc * (h-bfd))
        
        ! wetted perimiter, main channel, trapezoid
        WP = (bw + 2.0_prec * bfd * sqrt(1.0_prec + z*z))
        
        ! wetted perimiter, compound (floodplain), rectancular
        ! ??????????????????????????????????????????????????????????????????????????????????????????????????????????
        WPC = twcc + (2.0_prec*(h-bfd)) ! ??? I think this is wrong - should be (twcc - tw) + (2.0_prec*(h-bfd)) ???
        
        ! Hydraulic radius of main and compound channels combined
        R   = (AREA + AREAC)/(WP +WPC)  
        
    else ! if depth less than bankfull
        
        ! channel area, main channel, trapezoid
        AREA = (bw + h * z ) * h
        
        ! wetted perimiter , main channel, trapezoid
        WP = (bw + 2.0_prec * h * sqrt(1.0_prec + z*z))

        if(WP .gt. 0.0_prec) then
            R = AREA/WP
        else
            R = 0.0_prec
        endif
    endif

    !**kinematic celerity, c
    if(h .gt. bfd) then
    !*water outside of defined channel weight the celerity by the contributing area, and
    !*assume that the mannings of the spills is 2x the manning of the channel
        Ck = max(0.0_prec,((sqrt(s0)/n)*((5.0_prec/3.0_prec)*R**(2.0_prec/3.0_prec) - &
            ((2.0_prec/3.0_prec)*R**(5.0_prec/3.0_prec)*(2.0_prec*sqrt(1.0_prec + z*z)/(bw+2.0_prec*bfd*z))))*AREA &
            + ((sqrt(s0)/(ncc))*(5.0_prec/3.0_prec)*(h-bfd)**(2.0_prec/3.0_prec))*AREAC)/(AREA+AREAC))
    else
        if(h .gt. 0.0_prec) then !avoid divide by zero
            Ck = max(0.0_prec,(sqrt(s0)/n)*((5.0_prec/3.0_prec)*R**(2.0_prec/3.0_prec) - &
                ((2.0_prec/3.0_prec)*R**(5.0_prec/3.0_prec)*(2.0_prec*sqrt(1.0_prec + z*z)/(bw+2.0_prec*h*z)))))
        else
            Ck = 0.0_prec
        endif
    endif

    !**MC parameter, K
    
    ! K = dx/c, such Muskingum-Cunge method is equivalent to the finite diff form of the kinematic wave eqn.
    ! interesting if the travel time is shorter than the timestep, K == dt. Why?
    
    if(Ck .gt. 0.0_prec) then
        Km = max(dt,dx/Ck)
    else
        Km = dt
    endif

    !**MC parameter, X
    
    ! x = 1/2 - [D/(c*dx)], where D = Q/(2*tw*s0): Bedient & Huber p. 295
    ! interesting that X is set to 0.5 at most, and between 0.25 and 0.5 for the lower interval
    
    if(h .gt. bfd) then !water outside of defined channel
        if (interval .eq. upper_interval) then ! ???? upper interval? 
            X = min(0.5_prec,max(0.0_prec,0.5_prec*(1.0_prec-(Qj/(2.0_prec*twcc*s0*Ck*dx))))) 
        endif
        if (interval .eq. lower_interval) then ! ??? lower interval
            X = min(0.5_prec,max(0.25_prec,0.5_prec*(1.0_prec-(((C1*qup)+(C2*quc)+(C3*qdp) + C4)/(2.0_prec*twcc*s0*Ck*dx)))))
        endif
    else
        if(Ck .gt. 0.0_prec) then
            !H0
            if (interval .eq. upper_interval) then
                X = min(0.5_prec,max(0.0_prec,0.5_prec*(1.0_prec-(Qj/(2.0_prec*twl*s0*Ck*dx)))))
            endif
            !H
            if (interval .eq. lower_interval) then
                X = min(0.5_prec,max(0.25_prec,0.5_prec*(1.0_prec-(((C1*qup)+(C2*quc)+(C3*qdp) + C4)/(2.0_prec*twl*s0*Ck*dx)))))
            endif
        else
            X = 0.5_prec
        endif
    endif

    ! calculate MC D parameter 
    D = (Km*(1.0_prec - X) + dt/2.0_prec)              !--seconds
    if(D .eq. 0.0_prec) then
        !print *, "FATAL ERROR: D is 0 in MUSKINGCUNGE", Km, X, dt,D
        !call hydro_stop("In MUSKINGCUNGE() - D is 0.")
    endif

    ! calculate C1, C2, C3, C4 for the Cunge method: Bedient & Huber p. 295
    C1 =  (Km*X + dt/2.0_prec)/D
    C2 =  (dt/2.0_prec - Km*X)/D
    C3 =  (Km*(1.0_prec-X)-dt/2.0_prec)/D
    C4 =  (ql*dt)/D

    if (interval .eq. lower_interval) then
        if( (C4 .lt. 0.0_prec) .and. (abs(C4) .gt. (C1*qup)+(C2*quc)+(C3*qdp)))  then
            C4 = -((C1*qup)+(C2*quc)+(C3*qdp))
        endif
    endif
    !!Uncomment to show WP/WPC behavior above bankfull
    !if (interval .eq. upper_interval) then
    !    print *,"secant1 --", "WP:", WP, "WPC:", WPC
    !else
    !    print *,"secant2 --", "WP:", WP, "WPC:", WPC
    !endif

    if((WP+WPC) .gt. 0.0_prec) then  !avoid divide by zero
        Qj =  ((C1*qup)+(C2*quc)+(C3*qdp) + C4) - ((1.0_prec/(((WP*n)+(WPC*ncc))/(WP+WPC))) * &
                (AREA+AREAC) * (R**(2.0_prec/3.0_prec)) * sqrt(s0)) !f(x)
    else
        Qj = 0.0_prec
    endif

end subroutine secant2_h

!**---------------------------------------------------**!
!*                                                     *!
!*                 COURANT SUBROUTINE                  *!
!*                                                     *!
!**---------------------------------------------------**!
subroutine courant(h, bfd, bw, twcc, ncc, s0, n, z, dx, dt, ck, cn)

    implicit none

    real(prec), intent(in) :: h, bfd, bw, twcc, ncc, s0, n, z, dx, dt
    real(prec), intent(out) :: ck, cn
    
    real(prec) :: h_gt_bf, h_lt_bf, AREA, AREAC, WP, WPC, R
    
    h_gt_bf = max(h - bfd, 0.0_prec)
    h_lt_bf = min(bfd, h)
 
    AREA = (bw + h_lt_bf * z ) * h_lt_bf
    
    WP = (bw + 2 * h_lt_bf * sqrt(1 + z*z))
    
    AREAC = (twcc * h_gt_bf) 
    
    if(h_gt_bf .gt. 0.0_prec) then
        WPC = twcc + (2 * (h_gt_bf)) 
    else 
        WPC = 0
    endif
    
    R   = (AREA + AREAC)/(WP + WPC)
    
    ck = ((sqrt(s0)/n)* &
            ((5.0_prec/3.0_prec)*R**(2.0_prec/3.0_prec) - ((2.0_prec/3.0_prec)*R**(5.0_prec/3.0_prec)* &
            (2*sqrt(1.0_prec + z*z)/(bw+2.0_prec*h_lt_bf*z))))*AREA &
                + ((sqrt(s0)/(ncc))*(5.0_prec/3.0_prec)*(h_gt_bf)**(2.0_prec/3.0_prec))*AREAC)/(AREA+AREAC)
                
    
    cn = ck * (dt/dx)
  
end subroutine courant

end module muskingcunge_module