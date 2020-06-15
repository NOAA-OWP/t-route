module muskingcunge_module
    use var
    use precis
    implicit none

contains

    subroutine main
    !** This code implements Muskingum-Cunge method of NWM to a single reach that has
    !* multiple segments (nodes)
    !* refer to p.18-, WRF-HYDRO, RM1.
        integer :: i, j, k
        !*for times after initial time
        if (iseg>1) then
        !* take flow information from an upstream node of the same link(=reach)
        !* 1 indicates one previous time step
            qup=Qd(1,iseg-1)
            quc=Qd(1,iseg-1)
            qdp=Qd(1,iseg)
        elseif (uslinkflag==1) then
        !* At the first node of the current link j, when upstream link(=reach) is available,
        !* take flow information from the last node of the upstream link
            qup= Qus_prev
            quc= Qus_prev
            qdp= Qd(1,iseg)
        elseif (uslinkflag/=1) then
        !* At the first node of the current link j, when upstream link(=reach)
        !* is not available (when uslinkID<0), takes the following approach.
            qup=0.0_prec
            quc=0.0_prec
            qdp=Qd(1,iseg)
        end if
        qdc=0.0_prec
        vel=vela(1,iseg) !* 1 indicates one previous time step
        velp_chk = vel
        depth=deptha(1,iseg)
        depthp_chk = depth

        call mcNWM()

        Qd(2,iseg)= qdc !* 2 indicates current time step
        vela(2,iseg)= vel
        deptha(2,iseg)= depth

    end subroutine main
    !**---------------------------------------------------**!
    !*                                                     *!
    !*                  NEXT SUBROUTINE                    *!
    !*                                                     *!
    !**---------------------------------------------------**!
    subroutine mcNWM
        !* exactly follows SUBMUSKINGCUNGE in NWM:
        !* 1) qup and quc for a reach in upstream limit take zero values all the time
        !* 2) initial value of depth of time t of each reach is equal to the value at time t-1
        !* 3) qup as well as quc at time t for a downstream reach in a serial network takes
        !*    exactly the same value qdp at time t (or qdc at time t-1) for the upstream reach
        !use var_module2

        implicit none

        integer :: iter
        integer :: maxiter, tries
        real(prec) :: mindepth, aerror, rerror
        real(prec) :: R, Twl, h_1 , h, h_0, Qj, Qj_0

        !* parameters of Secant method
        maxiter  = 100
        mindepth=0.01_prec

        aerror = 0.01_prec
        rerror = 1.0_prec
        tries = 0

        if(Cs .eq. 0.00000000_prec) then
            z = 1.0_prec
        else
            z = 1.0_prec/Cs          !channel side distance (m)
        endif

        if(Bw .gt. Tw) then   !effectively infinite deep bankful
            bfd = Bw/0.00001_prec
        elseif (Bw .eq. Tw) then
            bfd =  Bw/(2.0_prec*z)  !bankfull depth is effectively
        else
            bfd =  (Tw - Bw)/(2.0_prec*z)  !bankfull depth (m)
        endif

        if (n .le. 0.0_prec .or. So .le. 0.0_prec .or. z .le. 0.0_prec .or. Bw .le. 0.0_prec) then
            !print*, "Error in channel coefficients -> Muskingum cunge", n, So, z, Bw
            !call hydro_stop("In MUSKINGCUNGE() - Error in channel coefficients")
        end if

        depth = max(depth, 0.0_prec)
        h     = (depth * 1.33_prec) + mindepth !1.50 of  depth
        h_0   = (depth * 0.67_prec)            !0.50 of depth

        if(ql .gt. 0.0_prec .or. qup .gt. 0.0_prec .or. qdp .gt. 0.0_prec) then  !only solve if there's water to flux
110 continue
            Qj_0 = 0.0_prec
            WPC = 0.0_prec
            AREAC = 0.0_prec
            iter = 0

            do while (rerror .gt. 0.01_prec .and. aerror .ge. mindepth .and. iter .le. maxiter)

                call secant_h0(h_0, Qj_0)

                call secant_h(h, Qj)

                if(Qj_0-Qj .ne. 0.0_prec) then
                    h_1 = h - ((Qj * (h_0 - h))/(Qj_0 - Qj)) !update h, 3rd estimate

                    if(h_1 .lt. 0.0_prec) then
                        h_1 = h
                    endif
                else
                    h_1 = h
                endif

                if(h .gt. 0.0_prec) then
                    rerror = abs((h_1 - h)/h) !relative error is new estatimate and 2nd estimate
                    aerror = abs(h_1 -h)      !absolute error
                else
                    rerror = 0.0_prec
                    aerror = 0.9_prec
                endif

                h_0  = max(0.0_prec,h)
                h    = max(0.0_prec,h_1)
                iter = iter + 1
                            !write(41,"(3i5,2x,8f15.4)") k, i, iter, dmy1, Qj_0, dmy2, Qj, h_0, h, rerror, aerror
                            !write(42,*) k, i, iter, dmy1, Qj_0, dmy2, Qj, h_0, h, rerror, aerror
                if( h .lt. mindepth) then  ! exit loop if depth is very small
                    goto 111
                endif
            end do !*do while (rerror .gt. 0.01 .and. ....
111    continue

            if(iter .ge. maxiter) then
                tries = tries + 1

                if(tries .le. 4) then  ! expand the search space
                    h     =  h * 1.33_prec
                    h_0   =  h_0 * 0.67_prec
                    maxiter = maxiter + 25 !and increase the number of allowable iterations
                    goto 110
                endif
                        !print*, "Musk Cunge WARNING: Failure to converge"
                        !print*, 'RouteLink index:', idx + linkls_s(my_id+1) - 1
                        !print*, "id,err,iters,tries",PC*nCC))/(WP+WPC))) * &
                        !        (AREA+AREAC) * (R**(2./3.)) * sqrt(So)) idx, rerror, iter, tries
                        !print*, "Ck,X,dt,Km",Ck,X,dt,Km
                        !print*, "So,dx,h",So,dx,h
                        !print*, "qup,quc,qdp,ql", qup,quc,qdp,ql
                        !print*, "bfd,Bw,Tw,Twl", bfd,Bw,Tw,Twl
                        !print*, "Qmc,Qmn", (C1*qup)+(C2*quc)+(C3*qdp) + C4,((1/(((WP*n)+(WPC*nCC))/(WP+WPC))) * &
                        !        (AREA+AREAC) * (R**(2./3.)) * sqrt(So))
            endif

            !*yw added for test
            !*DY and LKR Added to update for channel loss
            if(((C1*qup)+(C2*quc)+(C3*qdp) + C4) .lt. 0.0_prec) then
                if( (C4 .lt. 0.0_prec) .and. (abs(C4) .gt. (C1*qup)+(C2*quc)+(C3*qdp)) )  then ! channel loss greater than water in chan
                    qdc = 0.0_prec
                else
                    qdc = MAX( ( (C1*qup)+(C2*quc) + C4),((C1*qup)+(C3*qdp) + C4) )
                endif
            else
                qdc =  ((C1*qup)+(C2*quc)+(C3*qdp) + C4) !-- pg 295 Bedient huber
                !write(*,*)"C1", C1, "qup", qup, "C2", C2, "quc", quc, "C3", C3, "qdp", qdp, "C4", C4
            endif

            Twl = Bw + (2.0_prec*z*h)
            R = (h*(Bw + Twl) / 2.0_prec) / (Bw + 2.0_prec*(((Twl - Bw) / 2.0_prec)**2.0_prec + h**2.0_prec)**0.5_prec)
            vel =  (1.0_prec/n) * (R **(2.0_prec/3.0_prec)) * sqrt(So)  !*average velocity in m/s
            depth = h
        else   !*no flow to route
            qdc = 0.0_prec
            depth = 0.0_prec
        end if !*if(ql .gt. 0.0 .or. ...

    end subroutine mcNWM
    !**---------------------------------------------------**!
    !*                                                     *!
    !*                  NEXT SUBROUTINE                    *!
    !*                                                     *!
    !**---------------------------------------------------**!
    subroutine secant_h0(h_0, Qj_0)

        implicit none

        real(prec), intent(in) :: h_0  !, refQj_0
        real(prec), intent(out) :: Qj_0
        real(prec) :: Twl, AREA, WP, R, Ck, Km, X, D

        !**top surface water width of the channel inflow
        Twl = Bw + 2.0_prec*z*h_0

        !**hydraulic radius, R
        if(h_0 .gt. bfd) then !**water outside of defined channel
            AREA =  (Bw + bfd * z) * bfd
            AREAC = (TwCC * (h_0 -bfd)) !**assume compound component is rect. chan, that's 3 times the Tw
            WP = (Bw + 2.0_prec * bfd * sqrt(1.0_prec + z*z))
            WPC = TwCC + (2.0_prec * (h_0-bfd)) !**WPC is 2 times the Tw
            R   = (AREA + AREAC)/(WP +WPC)  !**hydraulic radius
            !print *, "warning: compound channel activated", h_0, bfd
        else
            AREA = (Bw + h_0 * z ) * h_0
            WP = (Bw + 2.0_prec * h_0 * sqrt(1.0_prec + z*z))
            !WPC = 0.0
            if(WP .gt. 0.0_prec) then
                R = AREA/ WP
            else
                R = 0.0_prec
            endif
        endif

        !**kinematic celerity, c
        if(h_0 .gt. bfd) then
        !*water outside of defined channel weight the celerity by the contributing area, and
        !*assume that the mannings of the spills is 2x the manning of the channel
            Ck = max(0.0_prec,((sqrt(So)/n)*((5.0_prec/3.0_prec)*R**(2.0_prec/3.0_prec) - &
                ((2.0_prec/3.0_prec)*R**(5.0_prec/3.0_prec)*(2.0_prec*sqrt(1.0_prec + z*z)/(Bw+2.0_prec*bfd*z))))*AREA &
                + ((sqrt(So)/(nCC))*(5.0_prec/3.0_prec)*(h_0-bfd)**(2.0_prec/3.0_prec))*AREAC)/(AREA+AREAC))
        else
            if(h_0 .gt. 0.0_prec) then
                Ck = max(0.0_prec,(sqrt(So)/n)*((5.0_prec/3.0_prec)*R**(2.0_prec/3.0_prec) - &
                    ((2.0_prec/3.0_prec)*R**(5.0_prec/3.0_prec)*(2.0_prec*sqrt(1.0_prec + z*z)/(Bw+2.0_prec*h_0*z)))))
            else
                Ck = 0.0_prec
            endif
        endif

        !**MC parameter, K
        if(Ck .gt. 0.000000_prec) then
            Km = max(dt,dx/Ck)
        else
            Km = dt
        endif

        !**MC parameter, X
        if(h_0 .gt. bfd) then !water outside of defined channel
            X = min(0.5_prec,max(0.0_prec,0.5_prec*(1.0_prec-(Qj_0/(2.0_prec*TwCC*So*Ck*dx)))))
            !X = min(0.5,max(0.0,0.5*(1-(refQj_0/(TwCC*So*Ck*dx)))))
        else
            if(Ck .gt. 0.0_prec) then
                X = min(0.5_prec,max(0.0_prec,0.5_prec*(1.0_prec-(Qj_0/(2.0_prec*Twl*So*Ck*dx)))))
                !X = min(0.5,max(0.0,0.5*(1-(refQj_0/(Twl*So*Ck*dx)))))
            else
                X = 0.5_prec
            endif
        endif

        !write(45,"(3i5,2x,4f10.3)") gk, gi, idx, h_0, Ck, Km, X
        D = (Km*(1.000_prec - X) + dt/2.0000_prec)              !--seconds
        if(D .eq. 0.0_prec) then
            !print *, "FATAL ERROR: D is 0 in MUSKINGCUNGE", Km, X, dt,D
            !call hydro_stop("In MUSKINGCUNGE() - D is 0.")
        endif

        C1 =  (Km*X + dt/2.000000_prec)/D
        C2 =  (dt/2.0000_prec - Km*X)/D
        C3 =  (Km*(1.00000000_prec-X)-dt/2.000000_prec)/D
        C4 =  (ql*dt)/D

        if((WP+WPC) .gt. 0.0_prec) then  !avoid divide by zero
            Qj_0 =  ((C1*qup)+(C2*quc)+(C3*qdp) + C4) - ((1.0_prec/(((WP*n)+(WPC*nCC))/(WP+WPC))) * &
                    (AREA+AREAC) * (R**(2.0_prec/3.0_prec)) * sqrt(So)) !f0(x)
        endif

    end subroutine secant_h0
    !**---------------------------------------------------**!
    !*                                                     *!
    !*                  NEXT SUBROUTINE                    *!
    !*                                                     *!
    !**---------------------------------------------------**!
    subroutine secant_h(h, Qj)

        implicit none

        real(prec),intent(in) :: h
        real(prec),intent(out) :: Qj
        real(prec) :: Twl, AREA, WP, R, Ck, Km, X, D

        !--upper interval -----------
        Twl = Bw + 2.0_prec*z*h  !--top width of the channel inflow

        if(h .gt. bfd) then !water outside of defined channel
            AREA =  (Bw + bfd * z) * bfd
            AREAC = (TwCC * (h-bfd)) !assume compound component is rect. chan, that's 3 times the Tw
            WP = (Bw + 2.0_prec * bfd * sqrt(1.0_prec + z*z))
            WPC = TwCC + (2.0_prec*(h-bfd)) !the additional wetted perimeter
            R   = (AREA + AREAC)/(WP +WPC)
            !print *, "warning: compound channel activated", h, bfd
        else
            AREA = (Bw + h * z ) * h
            WP = (Bw + 2.0_prec * h * sqrt(1.000000_prec + z*z))
            !WPC = 0.0
            if(WP .gt. 0.0_prec) then
                R = AREA/WP
            else
                R = 0.0_prec
            endif
        endif

        if(h .gt. bfd) then !water outside of defined channel, assumed rectangular, 3x TW and n = 3x
            Ck = max(0.0_prec,((sqrt(So)/n)*((5.0_prec/3.0_prec)*R**(2.0_prec/3.0_prec) - &
                ((2.0_prec/3.0_prec)*R**(5.0_prec/3.0_prec)*(2.0_prec*sqrt(1.0_prec + z*z)/(Bw + 2.0_prec*bfd*z))))*AREA &
                + ((sqrt(So)/(nCC))*(5.0_prec/3.0_prec)*(h-bfd)**(2.0_prec/3.0_prec))*AREAC)/(AREA+AREAC))
        else
            if(h .gt. 0.0_prec) then !avoid divide by zero
                Ck = max(0.0_prec,(sqrt(So)/n)*((5.0_prec/3.0_prec)*R**(2.0_prec/3.0_prec) - &
                    ((2.0_prec/3.0_prec)*R**(5.0_prec/3.0_prec)*(2.0_prec * sqrt(1.0_prec + z*z)/(Bw + 2.0_prec*h*z)))))
            else
                Ck = 0.0_prec
            endif
        endif

        if(Ck .gt. 0.0_prec) then
            Km = max(dt,dx/Ck)
        else
            Km = dt
        endif

        if(h .gt. bfd) then !water outside of defined channel
            X = min(0.5_prec,max(0.25_prec,0.5_prec*(1.0_prec-(((C1*qup)+(C2*quc)+(C3*qdp) + C4)/(2.0_prec*TwCC*So*Ck*dx)))))
        else
            if(Ck .gt. 0.0_prec) then
                X = min(0.5_prec,max(0.25_prec,0.5_prec*(1.0_prec-(((C1*qup)+(C2*quc)+(C3*qdp) + C4)/(2*Twl*So*Ck*dx)))))
            else
                X = 0.5_prec
            endif
        endif

        !write(45,"(3i5,2x,4f10.3)") gk, gi, idx, h, Ck, Km, X
        D = (Km*(1.0_prec - X) + dt/2.0_prec)              !--seconds
        if(D .eq. 0.0_prec) then
            !print *, "FATAL ERROR: D is 0 in MUSKINGCUNGE", Km, X, dt,D
            !call hydro_stop("In MUSKINGCUNGE() - D is 0.")
        endif

        C1 =  (Km*X + dt/2.000000_prec)/D
        C2 =  (dt/2.000000_prec - Km*X)/D
        C3 =  (Km*(1.000000_prec-X)-dt/2.000000_prec)/D
        C4 =  (ql*dt)/D

        if( (C4 .lt. 0.0_prec) .and. (abs(C4) .gt. (C1*qup)+(C2*quc)+(C3*qdp)))  then
            C4 = -((C1*qup)+(C2*quc)+(C3*qdp))
        endif

        if((WP+WPC) .gt. 0.0_prec) then
            Qj =  ((C1*qup)+(C2*quc)+(C3*qdp) + C4) - ((1.0000000_prec/(((WP*n)+(WPC*nCC))/(WP+WPC))) * &
                    (AREA+AREAC) * (R**(2.0_prec/3.0_prec)) * sqrt(So)) !f(x)
        endif

    end subroutine secant_h
end module muskingcunge_module

