module subtools

    !use constants
    !use arrays
    use var

contains

!+++-------------------------------------------------------------------------------
!+ computation of normal depth in regular/trapezoidal/compound channel using
!+ Newton-Raphson/Approximated/Bisection methods.
!+ Refer to Appendix C-2, Chaudhary and p71,RM1_MESH
!+++-------------------------------------------------------------------------------
    subroutine nmaldepcalc(ic, ynm0, dsc, ynm)

        implicit none

        integer, intent(in) :: ic
        real, intent(in) :: ynm0, dsc
        real, intent(out) :: ynm

        integer :: recnm_app, ib, ib_mx
        real :: So, sslp, bwd, manN
        real :: c0, c1
        real :: errNR, tolNR
        real :: lambda, epsr, zetn
        real :: yb0, yb1, yb2, fyn, yinc, ybp, mpy, cvrb, errb

        tolNR=0.01  !*tolerance
        c0=1.0    !*1 for SI units and 1.49 for customary English units
        lambda=c0
        manN= 1/sk(ic)
        !So= -(z(ic,jc) - z(ic-1,jc))/dx(ic-1,jc) !*this subroutine is called only in predictor step ** changed 5/13/2020
        So= chbtslp(ic-1)
        sslp= traps(ic)
        bwd=bo(ic)
        c1=manN*dsc/(c0*(So**0.5))
        tolNR=0.01  !*tolerance
        errNR=10.0*tolNR
        recnm_app=2     !* approach 1 or 2 for rectangular channel
        ynm=ynm0    !*initial estimate only for Newton-Rapson method.

        if (sslp==0.0) then
        !++-----------------------------
        !+ rectangular channel
        !++-----------------------------
            if (recnm_app==1) then
            !*Newton-Raphson method. Refer to Appendix C-2, Chaudhary and p71,RM1_MESH
            !*Determined not to be used.
            elseif (recnm_app==2) then
            !*Direct solution based on "Direct solutions for normal depth in parabolic and
            !*rectangular open channels using asymptotic matching technique"
                epsr=manN*dsc/(lambda*(bwd**(8.0/3.0))*(So**0.5))                         !*Eq.8
                zetn=((epsr**2.3733) + 4.3957*(epsr**3.0319) + 6.2205*(epsr**3.9556))**0.2528 !*Eq.26
                ynm=zetn*bwd                                                            !*Eq.7
            end if
        else
        !++--------------------------------------------------------------------------------------
        !+ normal depth either for overbank or inbank flow in trap.main ch. and rec. floodplains
        !+ using Bisection method, p95-96,RM1
        !++--------------------------------------------------------------------------------------
        !** Note: Even though the current depth (=used for ynm0) is of inbank flow, it is possible
        !*        to compute the final normal depth as in overbank flow. So, normal depth needs to
        !*        be computed considering both trapezoidal main and compound channels.
            fyn= fCCyn(ic, ynm0, dsc)
            if (fyn==0.0) then
                ynm=ynm0
                return
            elseif (fyn>0.0) then
                yb2=ynm0
                !* search for yb1
                yinc=0.1*ynm0
                ybp= ynm0 - yinc
                fyn= fCCyn(ic, ybp, dsc)
                do while ((fyn>0.0).and.(ybp>yinc))
                    ybp=ybp - yinc
                    fyn= fCCyn(ic, ybp, dsc)
                        !write(*,*) "yb1=",ybp, "fyn=",fyn
                enddo
                if (fyn==0.0) then
                    ynm= ybp
                    return
                elseif (fyn<0.0) then
                    yb1= ybp
                else
                    !* no available option then simply assume ynm is equal to the current depth.
                    ynm= ynm0
                    return
                endif
            else !* (fyn<0.0)
                yb1=ynm0
                !* search for yb2
                yinc=0.1*ynm0
                ybp= ynm0 + yinc
                fyn= fCCyn(ic, ybp, dsc)
                mpy=100.0  !*multiplier for y for deciding possibly largest value of y
                do while ((fyn<0.0).and.(ybp<mpy*ynm0))
                    ybp=ybp+yinc
                    fyn= fCCyn(ic, ybp, dsc)
                        !write(*,*) "yb2=",ybp, "fyn=",fyn
                enddo
                if (fyn==0.0) then
                    ynm= ybp
                    return
                elseif (fyn>0.0) then
                    yb2= ybp
                else
                    !* no available option then simply assume that ynm is equal to the current depth.
                    ynm= ynm0
                    return
                endif
            endif

            cvrb=0.01*abs((yb2-yb1)/2.0)
            ib_mx=20
            do ib=1,ib_mx
                yb0= (yb1 + yb2)/2.0
                fyn= fCCyn(ic, yb0, dsc)
                if (fyn==0.0) then
                    ynm= yb0
                    exit
                elseif (fyn<0.0) then
                    yb1=yb0
                else
                    yb2=yb0
                endif
                errb=abs(yb2-yb1)/2.0
                if (errb<cvrb) then
                    ynm=(yb1 + yb2)/2.0
                    exit
                end if
                    !write(*,*)"yb1=",yb1, "yb2=", yb2, "fyn=",fyn
            enddo
        endif

    end subroutine nmaldepcalc

!+++----------------------------------------------
!+ F(yn)=A**(5/3)/P**(2/3) - nQ/So**(1/2), p95,RM1
!+++----------------------------------------------
    real function fCCyn(ic, yxs, dsc)
        implicit none
        integer, intent(in) :: ic
        real, intent(in) :: yxs,dsc
        real :: bwd, sslp, Twd
        real :: ybf, TwCCi, A0,A1,A2, P0,P1,P2, tlA, tlP
        real :: So, manN, emanN, nom

        bwd= bo(ic)
        sslp= traps(ic)
        Twd= Tw(ic)
        !So= -(z(ic,jc) - z(ic-1,jc))/dx(ic-1,jc) !*this subroutine is called only in predictor step
        if (ic==1) then
            So= chbtslp(ic)        !** changed 5/13/2020
        else
            So= chbtslp(ic-1)        !** changed 5/13/2020
        endif

        ybf=(Twd - bwd)/(2.0*sslp)
        if (yxs<=ybf) then
        !* still inbank flow
            call areacalc(ic, yxs, A0)
            P0 = bwd + 2.0*yxs*(1.0 + sslp**2.0)**0.5
            manN= 1.0/sk(ic)
            fCCyn= A0**(5.0/3.0)/P0**(2.0/3.0) - manN*dsc/(So**0.5)
        else
        !* overbank flow
        !** A0,A1,A2 and P0,P1,P2 for overbank flow p84-2-1,RM1
            !* A0 & P0
            A0 = (bwd + sslp*ybf)*ybf + (yxs - ybf)*Twd
            P0 = bwd + 2.0*ybf*(1.0 + sslp**2.0)**0.5

            TwCCi=(TwCC(ic)-Twd)/2.0
            !* A1 & P1
            A1=(yxs - ybf)*TwCCi
            P1=TwCCi +  (yxs - ybf)
            !* A2 & P2
            A2= (yxs - ybf)*TwCCi
            P2=TwCCi +  (yxs - ybf)
            !** F(yn)
            tlA= A0+A1+A2
            tlP= P0+P1+P2
            !* equivalent Manning's N for compound channel, Eq(4-35),p91,Chaudhry
            nom= P0*(1.0/sk(ic))**(3.0/2.0) + P1*(1.0/skCC1(ic))**(3.0/2.0) + P2*(1.0/skCC2(ic))**(3.0/2.0)
            emanN= (nom/tlP)**(2.0/3.0)
            fCCyn= tlA**(5.0/3.0)/tlP**(2.0/3.0) - emanN*dsc/(So**0.5)
        endif

    end function fCCyn


!+++----------------------------------------------------------------------------
!+ computation of critical depth in regular/trapezoidal x-section using
!+ "Explicit solutions for critical and normal depths in trapezoidal and parabolic
!+  open channels" by Ali R.V.p69-70,RM1_MESH
!+++----------------------------------------------------------------------------
    !subroutine criticaldep_NR(sslp, bwd, dsc, ycr)
    subroutine criticaldep(ic, dsc, ycr)

        implicit none
        integer, intent(in) :: ic
        real, intent(in) :: dsc
        real, intent(out) :: ycr
        real :: sslp, bwd
        real :: dmy, epsc, tc0, etac

        sslp= traps(ic)
        bwd= bo(ic)
        if (sslp==0.0) then
        !* rectangular channel
            ycr=dsc**(2.0/3.0)/(grav*bwd**2)**(1.0/3.0)
        else
        !* trapezoidal channel
            !* approach 2 (based on "Explicit solutions for critical and normal depths
            !* in trapezoidal and parabolic open channels" by Ali R.V.
            epsc=4.0*((sslp**3.0)*(dsc**2.0)/(grav*bwd**5.0))**(1.0/3.0)

            dmy=(1.0+0.666*epsc**1.041)**0.374
            tc0=(1.0+1.161*epsc*dmy)**0.144

            dmy=(5.0*tc0**6.0 + 1)/(6.0*tc0**5.0 - epsc)
            !dmy=(5.0*tc0**0.864 + 1)/(6.0*tc0**0.72 - epsc)
            etac=-0.5+0.5*dmy**3.0
            ycr=bwd*etac/sslp
        endif
    end subroutine criticaldep

!+++----------------------------------------------------------------------------
!+ computation of downstream conjugate depth for regular/trapezoidal x-section
!+ when fow is supercritcal.
!+ Refer to p.40,Chaudhary and p72,RM1_MESH
!+++----------------------------------------------------------------------------
    subroutine conjugatedep(ic, ycj1, dsc, ycj2)
        implicit none

        integer, intent(in) :: ic
        real, intent(in) :: ycj1, dsc
        real, intent(out) :: ycj2
        integer iter, mxiter, trapcj_app
        real :: sslp, bwd, Twd, TCCwd
        real :: Axs, Fr1
        real :: uwQ, xcj, eta, kcj, delta, tcj, ycj0, lamX
        real :: thta, gma, pi1, mcj, lcj, zbar, ybf, nom, denm
        real :: acb, bcb, ccb, pcb, qcb, dscb, ycd1, ycd2, ycd3, arg1, arg2, pi314

        sslp= traps(ic)
        bwd= bo(ic)
        trapcj_app=2    !*approach 1 or 2 for trapezoidal channel

        if (sslp==0.0) then
        !* rectangular channel
            Axs=ycj1*bwd
            Fr1=dsc/((grav*(Axs**3.0)/bwd)**0.5)
            ycj2=0.5*ycj1*(-1.0 + (1.0+8.0*(Fr1**2.0))**0.5)
        else
        !** First, assume downstream conjugate depth is below bankfull depth
          !* inbank conjugate depth in trapezoidal channel
            !* Two approaches based on "Computing conjugate depths in trapezoidal channels" by Liu J. et al (2012).
            !* For more details, refer to p73~74, RM1_MESH
            uwQ=dsc/bwd                         !*unit width discharge
            xcj=ycj1/(uwQ**(2.0/3.0))
            eta=sslp*(uwQ**(2.0/3.0))/bwd         !*eta
            kcj=4.0*eta*((1.0/grav)**(1.0/3.0))

            delta=(((1.0+kcj*((1.0+kcj)**0.2))**0.5)-1.0)/2.0
            tcj=delta/eta

            ycj0=tcj + (tcj-xcj)*((tcj/xcj)**0.5)
            lamX=6.0/(grav*xcj*(1+eta*xcj))+(xcj**2.0)*(3.0+2.0*eta*xcj)

            if (trapcj_app==1) then
            !* Approach 1 (iterative formula)
                mxiter=3
                do iter=1,mxiter
                    ycj0=((lamX-6.0/(grav*ycj0*(1.0+eta*ycj0)))/(3.0+2.0*eta*ycj0))**0.5
                end do
            elseif (trapcj_app==2) then
            !* Approach 2 (estimation formula)
                ycj0=((lamX-6.0/(grav*ycj0*(1.0+eta*ycj0)))/(3.0+2.0*eta*ycj0))**0.5
            end if
            ycj2=ycj0*(uwQ**(2.0/3.0))

            !** When downstream conjugate depth turns out to be higher than bankfull depth,
            !** apply compound channel solution.
            Twd= Tw(ic)
            ybf= (Twd - bwd)/(2.0*sslp)
            if (ycj2>ybf) then
            !* compound channel with trapezoidal main and rec. floodplains
                !Twd= Tw(ic,jc)
                TCCwd= TwCC(ic)
                !* conjugate depth based Eq(2.54),p40,Chaudhry. Also, refer to p94,RM1
                ybf= (Twd - bwd)/(2.0*sslp)
                mcj= ybf*(3.0*bwd + 4.0*sslp*ybf)/(6.0*(bwd + sslp*ybf))
                lcj= ybf*(bwd + sslp*ybf)
                !* zbar at location 1 before hydraulic jump, p94,RM1
                nom= 0.5*TCCwd*(ycj1**2.0 + ybf**2.0) + (lcj- TCCwd*ybf)*ycj1 - mcj*lcj
                denm= TCCwd*(ycj1 - ybf) + lcj
                zbar= nom/denm
                call areacalc(ic, ycj1, Axs)
                !* pi at location 1, p94-1,RM1
                pi1= dsc**2.0/(grav*Axs) + zbar*Axs
                gma= 0.5*TCCwd*ybf**2.0 - mcj*lcj
                thta= lcj - TCCwd*ybf
                !* parameters of cubic poly. p94-1,RM1
                acb= 3.0*thta/TCCwd
                bcb= 2.0*(thta**2.0/TCCwd + gma - pi1)/TCCwd
                ccb= 2.0*(thta*(gma-pi1) + dsc**2.0/grav)/TCCwd**2.0
                !* Using solution in "Solving cubic equations". Also, p94-2~-3,RM1
                pcb= bcb - acb**2.0/3.0
                qcb= 2.0*acb**3.0/27.0 - acb*bcb/3.0 + ccb
                dscb= qcb**2.0/4.0 + pcb**3.0/27.0

                if (dscb>0.0) then
                !* only one real soln.
                    ycj2= (-qcb/2.0 + dscb**0.5)**(1.0/3.0) + (-qcb/2.0 - dscb**0.5)**(1.0/3.0) - acb/3.0
                elseif (dscb==0.0) then
                    ycd1= -2.0*(qcb/2.0)**(1.0/3.0) - acb/3.0
                    ycd2= (qcb/2.0)**(1.0/3.0) - acb/3.0
                    ycj2= max(ycd1, ycd2, ycj1)
                else
                    arg1= 3.0*3.0**0.5*qcb/(2.0*((-pcb)**0.5)**3.0)
                    !* inverse of sine func approximated by infinite series, p94-2,RM1
                    arg2= arg1 + arg1**3.0/6.0 + 3.0*arg1**5.0/40.0 + 5.0*arg1**7.0/112.0
                    ycd1= (2.0/3.0**0.5)*(-pcb)**0.5*sin(arg2/3.0) - acb/3.0
                    pi314= 3.141593
                    ycd2= -(2.0/3.0**0.5)*(-pcb)**0.5*sin(arg2/3.0 + pi314/3.0) - acb/3.0
                    ycd3= (2.0/3.0**0.5)*(-pcb)**0.5*sin(arg2/3.0 + pi314/6.0) - acb/3.0
                    ycj2= max(ycd1, ycd2, ycd3, ycj1)
                endif
            endif
        endif

    end subroutine conjugatedep

!+++-----------------------------------------------------
!+ Computation of depth with given area, p77, RM1_MESH
!+++-----------------------------------------------------
    subroutine depthcalc(ic, Axs, yxs)
        implicit none

        integer, intent(in) :: ic
        real, intent(in) :: Axs
        real, intent(out) :: yxs
        real :: bwd, sslp, ybf, Abf

        bwd= bo(ic)
        sslp= traps(ic)

        if (sslp==0.0) then
        !* rectangular channel
            yxs=Axs/bo(ic)
        else
            ybf= (Tw(ic) - bwd)/(2.0*sslp)
            !call areacalc(ic, ybf, Abf)
            Abf= (bwd + sslp*ybf)*ybf
            if (Axs<=Abf) then
            !* inbank flow in trapezoidal channel
                yxs=(- bwd + (bwd**2.0 + 4.0*sslp*Axs)**0.5)/(2.0*sslp)
            else
            !* overbank on compound channel having trapezoidal main and rectangular floodplains, p80-,RM1
                yxs=ybf + (Axs - Abf)/TwCC(ic)
            end if
        endif
    end subroutine depthcalc

!+++-----------------------------------------------------
!+ Computation of depth with given area and channel data
!+++-----------------------------------------------------
    subroutine depthcalc2(bxs, ssxs, twxs, twccxs, Abf, Axs, yxs)
    !subroutine depthcalc2(bxs, ssxs, twxs, twccxs, chbtslpxs, Abf, Axs, yxs)
        implicit none
        real, intent(in) :: bxs, ssxs, twxs, twccxs, Abf, Axs
        real, intent(out) :: yxs
        real :: ybf

        if (ssxs==0.0) then
        !* rectangular channel
            yxs=Axs/bxs
        else
            if (Axs<=Abf) then
            !* inbank flow in trapezoidal channel
                yxs=(- bxs + (bxs**2.0 + 4.0*ssxs*Axs)**0.5)/(2.0*ssxs)
            else
            !* overbank on compound channel having trapezoidal main and rectangular floodplains, p80-,RM1
                ybf= (twxs - bxs)/(2.0*ssxs)
                yxs= ybf + (Axs - (bxs + ssxs*ybf)*ybf)/twccxs
            end if
        endif

    end subroutine depthcalc2

!+++-----------------------------------------------------
!+ Computation of depth with given area and channel data
!+++-----------------------------------------------------
    subroutine depthcalc3(bxs, ssxs, twxs, twccxs, Axs, yxs)

        implicit none
        real, intent(in) :: bxs, ssxs, twxs, twccxs, Axs
        real, intent(out) :: yxs
        real :: ybf, Abf !, yxs

        if (ssxs==0.0) then
        !* rectangular channel
            yxs=Axs/bxs
        else
            ybf= (twxs - bxs)/(2.0*ssxs)
            Abf= (bxs + ssxs*ybf)*ybf
            if (Axs<=Abf) then
            !* inbank flow in trapezoidal channel
                yxs=(- bxs + (bxs**2.0 + 4.0*ssxs*Axs)**0.5)/(2.0*ssxs)
            else
            !* overbank on compound channel having trapezoidal main and rectangular floodplains, p80-,RM1
                ybf= (twxs - bxs)/(2.0*ssxs)
                yxs= ybf + (Axs - (bxs + ssxs*ybf)*ybf)/twccxs
            end if
        endif
        !sbtoolrst= yxs !*make this result accessible globally.

    end subroutine depthcalc3

!+++---------------------------------------------------------------------------
!+ Computation of area of various channel x-sections with depth as an argument
!+ and without pre-specified chshp
!+++---------------------------------------------------------------------------
    subroutine areacalc(ic, yxs, Axs)

        implicit none

        integer, intent(in) :: ic
        real, intent(in) :: yxs
        real, intent(out) :: Axs
        real :: bwd, sslp, ybf

        bwd= bo(ic)
        sslp= traps(ic)

        if (sslp==0.0) then
        !* rectangular channel
            Axs=bwd*yxs
        else
        !* determine whether inbank or overbank flow
            !* bankfull depth
            ybf= (Tw(ic) - bwd)/(2.0*sslp)
            if (yxs.le.ybf) then
                !* trapezoidal channel inbank flow
                Axs=(bwd + sslp*yxs)*yxs
            else
                !*overbank flow on rect. floodplains
                Axs=(bwd + sslp*ybf)*ybf + TwCC(ic)*(yxs-ybf)
            end if
        endif

    end subroutine areacalc

!+++---------------------------------------------------------------------------
!+ Computation of area of various channel x-sections with depth as an argument
!+ and without pre-specified chshp
!+++---------------------------------------------------------------------------
    subroutine areacalc2(bxs, ssxs, twxs, twccxs, yxs, Axs)

        implicit none

        !integer, intent(in) :: ic
        real, intent(in) :: bxs, ssxs, twxs, twccxs, yxs
        real, intent(out) :: Axs
        real :: ybf !bwd, sslp,

        !bwd= bo(ic)
        !sslp= traps(ic)
        if (ssxs==0.0) then
        !* rectangular channel
            Axs=bxs*yxs
        else
        !* determine whether inbank or overbank flow
            !* bankfull depth
            ybf= (twxs - bxs)/(2.0*ssxs)
            if (yxs.le.ybf) then
                !* trapezoidal channel inbank flow
                Axs=(bxs + ssxs*yxs)*yxs
            else
                !*overbank flow on rect. floodplains
                Axs=(bxs + ssxs*ybf)*ybf + twccxs*(yxs-ybf)
            end if
        endif
        !sbtoolrst= Axs !* to make Axs globally accessible
    end subroutine areacalc2

!+++-----------------------------------------------------
!+ Computation of I1
!+++-----------------------------------------------------
    subroutine I1calc(ic, dp, Axs, I1)

        implicit none

        integer, intent(in) :: ic
        real, intent(in) :: dp, Axs
        real, intent(out) :: I1
        real :: I1m, I1f, ybf, Abf
        real :: sslp, bwd

        bwd= bo(ic)
        sslp= traps(ic)

        if (sslp==0.0) then
        !* rectangular channel
            !*  p.80-1, RM1
            I1=0.5*bwd*(dp**2.0)
        else
            ybf= (Tw(ic) - bwd)/(2.0*sslp)
            call areacalc(ic, ybf, Abf)
            if (Axs<=Abf) then
            !* inbank flow in trapezoidal channel
                !* p.80-3, RM1
                I1=(dp**2.0)*(3.0*bwd + 2.0*sslp*dp)/6.0
            else
            !* overbank flow on compound channel having trapezoidal main and rectangular floodplains, p80-5,RM1
                I1m=ybf*(bwd*dp + (sslp*dp-0.5*bwd)*ybf - (2.0*sslp/3.0)*ybf**2.0)
                I1f=0.5*TwCC(ic)*(dp-ybf)**2.0
                I1=I1m + I1f
            endif
        endif

    end subroutine I1calc


!+++-----------------------------------------------------
!+ Computation of hydraulic radius R (=A/P)
!+++-----------------------------------------------------
    subroutine hydRcalc(ic, dp, Axs, hydR)

        implicit none

        integer, intent(in) :: ic
        real, intent(in) :: dp, Axs
        real, intent(out) :: hydR
        real :: bwd, sslp, Twd, ybf, Abf

        bwd= bo(ic)
        sslp= traps(ic)
        Twd= Tw(ic)

        if (sslp==0.0) then
        !* rectangular channel
            hydR=Axs/(bwd + 2.0*dp)
        else
            ybf= (Twd - bwd)/(2.0*sslp)
            call areacalc(ic, ybf, Abf)
            if (Axs<=Abf) then
            !* inbank flow in trapezoidal channel
                hydR=Axs/(bwd + 2.0*dp*((1.0+(sslp**2.0))**0.5))
            else
            !* compound channel having trapezoidal main and rectangular floodplains, p80-,RM1
                hydR=Axs/(bwd + 2.0*ybf*(1.0 + sslp**2.0)**0.5 + TwCC(ic) - Twd + 2.0*(dp-ybf))
            endif
        endif

    end subroutine hydRcalc

!+++-----------------------------------------------------
!+ Computation of hydraulic radius R (=A/P)
!+++-----------------------------------------------------
    subroutine hydRcalc2(bxs, ssxs, twxs, twccxs, yxs, Abf, Axs, hydR)
        implicit none
        real, intent(in) :: bxs, ssxs, twxs, twccxs, yxs, Abf, Axs
        real, intent(out) :: hydR
        real :: ybf

        if (ssxs==0.0) then
        !* rectangular channel
            !hydR=Axs/(bxs + 2.0*yxs)
            hydR= Axs/(bxs + 2.0*yxs)
        else
            !ybf= (Twd - bwd)/(2.0*sslp)
            ybf= (twxs - bxs)/(2.0*ssxs)
            !call areacalc(ic, jc, ybf, Abf)
            if (Axs<=Abf) then
            !* inbank flow in trapezoidal channel
                !hydR=Axs/(bwd + 2.0*dp*((1.0+(sslp**2.0))**0.5))
                hydR=Axs/(bxs + 2.0*yxs*((1.0 + ssxs**2.0)**0.5))
            else
            !* compound channel having trapezoidal main and rectangular floodplains, p80-,RM1
                !hydR=Axs/(bwd + 2.0*ybf*(1.0 + sslp**2.0)**0.5 + TwCC(ic,jc) - Twd + 2.0*(dp-ybf))
                hydR= Axs/(bxs + 2.0*ybf*(1.0 + ssxs**2.0)**0.5 + twccxs - twxs + 2.0*(yxs-ybf))
            endif
        endif

    end subroutine hydRcalc2


!+++-----------------------------------------------------
!+ Computation of conveyance K (=1/N * A * R^(2/3))
!+++-----------------------------------------------------
    subroutine Kcalc(ic, dp, Axs, hydR, cnvey)

        implicit none

        integer, intent(in) :: ic
        real, intent(in) :: dp, Axs, hydR
        real, intent(out) :: cnvey
        real :: bwd, sslp, Twd, TCCwd
        real :: subA, subP, ybf, TwCCi, K0, K1, K2

        bwd= bo(ic)
        sslp= traps(ic)
        !* top width of main channel
        Twd=Tw(ic)
        !* top width of compound channel
        TCCwd=TwCC(ic)
        !* bankfull depth
        ybf= (Twd - bwd)/(2.0*sslp)

        if ((sslp==0.0).or.(dp<=ybf)) then
        !* inbank flow in rectangular or trapezoidal channel
            cnvey=sk(ic)*Axs*(hydR**(2.0/3.0))
        else
        !* overbank flow in compound channel having trapezoidal main channel with rectangular floodplains, p84-2~p84-2-1, RM1
            !* conveyance in the main channel, K0
            subA = (bwd + sslp*ybf)*ybf + (dp - ybf)*Twd
            subP = bwd + 2.0*ybf*(1.0+ sslp**2.0)**0.5
            K0 = sk(ic)*subA**(5.0/3.0)/subP**(2.0/3.0)
            !* conveyance in the left floodplain (assuming symmetric), K1
            TwCCi=(TCCwd - Twd)/2.0
            subA=(dp - ybf)*TwCCi
            subP=TwCCi +  (dp - ybf)
            K1 = skCC1(ic)*subA**(5.0/3.0)/subP**(2.0/3.0)
            !* conveyance in the left floodplain (assuming symmetric), K1
            TwCCi=(TCCwd - Twd)/2.0
            subA=(dp - ybf)*TwCCi
            subP=TwCCi +  (dp - ybf)
            K2 = skCC2(ic)*subA**(5.0/3.0)/subP**(2.0/3.0)

            cnvey=K0+K1+K2
        endif

    end subroutine Kcalc

!+++-------------------------------------------------------------------------------------
!+ Computation of db/dx for inbank flow at current node, p81-2,RM1_MESH
!+++-------------------------------------------------------------------------------------
    subroutine dbdxcalc(ic, iup, idw, dp, dbdxval)
        implicit none

        integer, intent(in) :: ic, iup, idw
        real, intent(in) :: dp
        real, intent(out) :: dbdxval
        real :: sslp

        sslp= traps(ic)
        if (sslp==0.0) then
        !* For inbank flow in rectangular (traps()=0.0)
            dbdxval=(bo(idw)-bo(iup))/dx(iup)
        else
        !* For inbank flow in trapezoidal channel
            dbdxval= (bo(idw) + 2.0*traps(idw)*dp - (bo(iup)+ 2.0*traps(iup)*dp))/dx(iup)
        end if

    end subroutine dbdxcalc

!+++-----------------------------------------------------
!+ Computation of dbdx for overbank flow at current node
!+++-----------------------------------------------------
    subroutine dbdxCCcalc(ic, iup, idw, dbdxm, dbdxf)
        implicit none

        integer, intent(in) :: ic, iup, idw
        real, intent(out) :: dbdxm, dbdxf
        real :: ybf

        ybf= (Tw(ic) - bo(ic))/(2.0*traps(ic))
        !* For overbank flow at current node
        dbdxm= (bo(idw) + 2.0*traps(idw)*ybf - (bo(iup)+ 2.0*traps(iup)*ybf))/dx(iup)
        dbdxf= (TwCC(idw)-TwCC(iup))/dx(iup)

    end subroutine dbdxCCcalc

!+++-----------------------------------------------------
!+ I2 computation for compound channel at current node
!+++-----------------------------------------------------
    subroutine I2CCcalc(ic, dbdxm, dbdxf, dp, I2)
        implicit none

        integer, intent(in) :: ic
        real, intent(in) :: dbdxm, dbdxf, dp
        real, intent(out) :: I2
        real :: ybf

        !* compound channel having trapezoidal main channel with rectangular floodplain, p81-, p42-c, p43-c,RM1
        !* bankfull depth
        ybf=(Tw(ic) - bo(ic))/(2.0*traps(ic))
        !* For overbank flow
        I2= dbdxm*ybf*(dp - 0.5*ybf) + 0.5*dbdxf*(dp - ybf)**2.0

    end subroutine I2CCcalc

!+++-----------------------------------------------------
!+ I2 computation, p.114~115, RM3
!+++-----------------------------------------------------
    subroutine I2calc(ic, iup, idw, delx, hc, I2)
        implicit none
        integer, intent(in) :: ic, iup, idw
        real, intent(in) :: delx, hc
        real, intent(out) :: I2
        real :: mb, ms, mtwcc
        real :: hcbf

        !* change of bottom width in x
        mb= (bo(idw) - bo(iup))/delx
        !* change of side-slope in x
        ms= (traps(idw) - traps(iup))/delx
        !* change of top width of floodplain in x
        mtwcc= (TwCC(idw) - TwCC(iup))/delx

        if (traps(ic)==0.0) then
            !* rect. main ch. flow
            I2= mb
        else
            hcbf=(Tw(ic) - bo(ic))/(2.0*traps(ic))
            if (hc.le.hcbf) then
            !* trapezoidal main ch. flow, p114,RM3
                I2= (0.5*mb + ms*hc/3.0)*hc**2.0
            else
            !* compound ch. flow, p115,RM3
                I2= 0.5*mb*hcbf*(2.0*hc - hcbf) + 2.0*ms*(hcbf**2.0)*(3.0*hc/2.0-hcbf)/3.0 + 0.5*mtwcc*(hc-hcbf)**2.0
            endif
        endif

    end subroutine I2calc

!+++-----------------------------------------------------
!+ I2 computation, p.103-3~-5, RM4
!+++-----------------------------------------------------
!    subroutine I2calc(ic, iup, idw, delx, ho, hc, I2)
!
!        implicit none
!        integer, intent(in) :: ic, iup, idw, jc
!        real, intent(in) :: delx, ho, hc
!        real, intent(out) :: I2
!        real :: mtw, mtwcc
!        real :: hbfc, hbfup, hbfdw
!
!        !* change of top width of trap.main in x
!        mtw= bo(idw) + 2.0*traps(idw)*ho  - (bo(iup) + 2.0*traps(iup)*ho)
!        mtw= mtw/delx
!        !* change of top width of floodplain in x
!        mtwcc= (TwCC(idw) - TwCC(iup))/delx
!
!        hbfup=(Tw(iup) - bo(iup))/(2.0*traps(iup))
!        hbfdw=(Tw(idw) - bo(idw))/(2.0*traps(idw))
!
!        if ((ho<=hbfup).or.(ho<=hbfdw)) then
!        !* rec/trapezoidal main ch. flow, p.103-3, RM4
!            I2= 0.5*mtw*(hc**2.0)
!        else
!        !* compound ch. flow, p.103-4, RM4
!            I2= 0.5*mtwcc*(hc**2.0)
!        endif
!
!
!    end subroutine I2calc

!+++-----------------------------------------------------
!+ Computation of c in matrix A
!+++-----------------------------------------------------
    subroutine c_mtrxAcalc(ic, yxs, Axs, cA)
        implicit none

        integer, intent(in) :: ic
        real, intent(in) :: yxs,  Axs
        real, intent(out) :: cA
        real :: bwd, sslp, nom, denom, ybf, term1, term2

        bwd= bo(ic)
        sslp= traps(ic)

        if (sslp==0.0) then
        !* inbank flow in rectangular channel
            !* p80-2, RM1
            cA=(grav*Axs/bwd)**0.5
        else
            ybf=(Tw(ic) - bwd)/(2.0*sslp)
            if (yxs<=ybf) then
           !* inbank flow in trapezoidal channel, p80-4-2,RM1
                nom=yxs*(bwd + sslp*yxs)
                denom=(bwd**2.0 + 4.0*sslp*Axs)**0.5
                cA=(grav*nom/denom)**0.5
            else
           !* overbank flow in compound ch. with trap. main and rect. floodplains, p80-6,RM1
                term1=grav*ybf*(bwd + sslp*ybf)/TwCC(ic)
                term2=grav*(yxs-ybf)
                cA=(term1+term2)**0.5
            endif
        endif

    end subroutine c_mtrxAcalc

!+++----------------------------------------------------------
!+ Computation of dI1/dA for c in matrix A based on p112~3, RM3
!+++----------------------------------------------------------
    subroutine dI1dAcalc(ic, hc, dI1dA)
    !* dI1dA calculates based on equations in p.112~3,RM3
        implicit none
        integer, intent(in) :: ic
        real, intent(in) :: hc
        real, intent(out) :: dI1dA
        real :: bwd, sslp, twd, tccwd
        real :: hcbf

        if (traps(ic)==0.0) then
            !* rect. main ch. flow
            dI1dA= hc
        else
            bwd= bo(ic)
            sslp= traps(ic)
            twd= Tw(ic)
            tccwd= TwCC(ic)
            hcbf=(twd - bwd)/(2.0*sslp)

            if (hc.le.hcbf) then
            !* trapezoidal main ch. flow, p112,RM3
                dI1dA= (bwd + sslp*hc)*hc/(bwd + 2.0*sslp*hc)
            else
            !* compound ch. flow, p113,RM3
               dI1dA= ((bwd + sslp*hcbf - tccwd)*hcbf + tccwd*hc)/tccwd
            endif
        endif

    end subroutine dI1dAcalc

!+++------------------------------------------------------------------------
!+ Computation of dI2/dA for (2,1) component of matrix S based on p114~5, RM3
!+++------------------------------------------------------------------------
     subroutine dI2dAcalc(ic, iup, idw, delx, hc, dI2dA)
    !* dI2dA calculates based on equations in p.I1-1~9,RM3
        implicit none
        integer, intent(in) :: ic, iup, idw
        real, intent(in) :: delx, hc
        real, intent(out) :: dI2dA
        real :: bwd, sslp, twd, tccwd, hcbf
        real :: bwdup, bwddw, mb, sslpup, sslpdw, ms, tccwdup, tccwddw, mtwcc

        bwd= bo(ic)
        sslp= traps(ic)
        twd= Tw(ic)
        tccwd= TwCC(ic)
        !* change of bottom width in x
        bwdup= bo(iup)
        bwddw= bo(idw)
        mb= (bwddw - bwdup)/delx
        !* change of side-slope in x
        sslpup= traps(iup)
        sslpdw= traps(idw)
        ms= (sslpdw - sslpup)/delx
        !* change of top width of floodplain in x
        tccwdup= TwCC(iup)
        tccwddw= TwCC(idw)
        mtwcc= (tccwddw - tccwdup)/delx

        if (traps(ic)==0.0) then
            !* rect. main ch. flow
            dI2dA= mb*hc/bwd
        else
            hcbf=(twd - bwd)/(2.0*sslp)
            if (hc.le.hcbf) then
            !* trapezoidal main ch. flow, p114,RM3
                dI2dA= (mb + ms*hc)*hc/(bwd + 2.0*sslp*hc)
            else
            !* compound ch. flow, p115,RM3
               dI2dA= ((mb + ms*hcbf)*hcbf + mtwcc*(hc-hcbf))/tccwd
            endif
        endif
   end subroutine dI2dAcalc

    !* chmay
!+++------------------------------------------------------------------------
!+ Computation of dK/dA for (2,1) component of matrix S based on 116~8, RM3
!+++------------------------------------------------------------------------
     subroutine dKdAcalc_p116(ic, hc, dKdA)

        implicit none
        integer, intent(in) :: ic
        real, intent(in) :: hc
        real, intent(out) :: dKdA
        real :: bwd, sslp, twd, tccwd, hcbf
        real :: Ac, Rc, manN, dkda1, P0, P1, P2, tlP,nom, emanN

        bwd= bo(ic)
        sslp= traps(ic)
        twd= Tw(ic)
        tccwd= TwCC(ic)

        if (traps(ic)==0.0) then
            !* rect. main ch. flow
            call areacalc(ic, hc, Ac)
            call hydRcalc(ic, hc, Ac, Rc)
            manN=1.0/sk(ic)

            dKdA= (Rc**(2.0/3.0))*(1.0 + 2.0*(1.0 - 2.0*Rc/bwd)/3.0)/manN
        else
            hcbf=(twd - bwd)/(2.0*sslp)
            if (hc.le.hcbf) then
            !* trapezoidal main ch. flow, p117,RM3
                call areacalc(ic, hc, Ac)
                call hydRcalc(ic, hc, Ac, Rc)
                manN=1.0/sk(ic)

                dkda1= 2.0*Rc*((1.0+sslp**2.0)**0.5)/(bwd + 2.0*sslp*hc)
                dKdA= (Rc**(2.0/3.0))*(1.0 + (2.0/3.0)*(1.0 - dkda1))/manN
            else
                !* compound ch. flow, p118,RM3
                call areacalc(ic, hc, Ac)
                call hydRcalc(ic, hc, Ac, Rc)
                !* equivalent Manning's N for compound channel, Eq(4-35),p91,Chaudhry
                P0 = bwd + 2.0*hcbf*(1.0 + sslp**2.0)**0.5
                P1= (tccwd-twd)/2.0 +  (hc - hcbf)
                P2= (tccwd-twd)/2.0 +  (hc - hcbf)
                tlP= P0+P1+P2
                nom= P0*(1.0/sk(ic))**(3.0/2.0) + P1*(1.0/skCC1(ic))**(3.0/2.0) + P2*(1.0/skCC2(ic))**(3.0/2.0)
                emanN= (nom/tlP)**(2.0/3.0)

                dKdA= (Rc**(2.0/3.0))*(1.0 + 2.0*(1.0 - 2.0*Rc/tccwd)/3.0)/emanN
            endif
        endif
     end subroutine dKdAcalc_p116

!+++-----------------------------------------------------
!+ Computation of dk/dA,p82-83,RM1_MESH
!+++-----------------------------------------------------
    !subroutine dKdAcalc(manN, sslp, bwd, Axs, dkda)
    subroutine dKdAcalc(ic, yxs, Axs, dkda)
        implicit none
        integer, intent(in) :: ic
        real, intent(in) ::  yxs, Axs
        real, intent(out) :: dkda
        real :: eta, pi, phi, bwd, sslp, manN, Twd, TCCwd
        real :: ybf, wpi, Ai, TwCCi, trm1, trm2, dkda0, dkda1, dkda2

        bwd= bo(ic)
        sslp= traps(ic)
        manN=1.0/sk(ic)

        if (sslp==0.0) then
        !* rectangular channel
            eta=(bwd + 2.0*Axs/bwd)
            dkda=((5.0/3.0)*(Axs**(2.0/3.0))*eta - 4.0*Axs**(5.0/3.0)/(3.0*bwd))/(manN*(eta**(5.0/3.0)))
        else
            ybf=(Tw(ic) - bwd)/(2.0*sslp)
            if (yxs<=ybf) then
            !* inbank flow in trapezoidal channel
                pi= bwd + (-bwd + (bwd**2.0 + 4.0*sslp*Axs)**0.5)*&
                    ((1.0+(sslp**2.0))**0.5)/sslp
                phi=2.0*((1.0 + sslp**2.0)**0.5)*((bwd**2.0 + 4.0*sslp*Axs)**(-0.5))
                dkda=((5.0/3.0)*(Axs**(2.0/3.0))*pi - (2.0/3.0)*(Axs**(5.0/3.0))*phi)/(manN*(pi**(5.0/3.0)))
            else
            !* overbank flow in compound ch. with trap. main and rect. floodplains, p84-3,RM1
                !* dKo/dA
                Twd= Tw(ic)
                TCCwd= TwCC(ic)
                ybf= (Twd - bwd)/(2.0*sslp)
                wpi= bwd + 2.0*ybf*(1.0 + sslp**2.0)**0.5
                Ai= (bwd + sslp*ybf)*ybf + (yxs-ybf)*Twd
                trm1= 5.0*Twd*Ai**(2.0/3.0)/(3.0*TCCwd*wpi**(2.0/3.0))
                trm2= 4.0*Ai**(5.0/3.0)/(3.0*TCCwd*wpi**(5.0/3.0))
                dkda0=(trm1 - trm2)/manN
                !* dK1/dA
                manN=1.0/skCC1(ic)
                TwCCi=(TCCwd - Twd)/2.0
                wpi= TwCCi + (yxs-ybf)
                Ai= (yxs-ybf)*TwCCi
                trm1= 5.0*TwCCi*Ai**(2.0/3.0)/(3.0*TCCwd*wpi**(2.0/3.0))
                trm2= 4.0*Ai**(5.0/3.0)/(3.0*TCCwd*wpi**(5.0/3.0))
                dkda1= (trm1 - trm2)/manN
                !* dK2/dA
                manN=1.0/skCC2(ic)
                trm1= 5.0*TwCCi*Ai**(2.0/3.0)/(3.0*TCCwd*wpi**(2.0/3.0))
                trm2= 4.0*Ai**(5.0/3.0)/(3.0*TCCwd*wpi**(5.0/3.0))
                dkda2= (trm1 - trm2)/manN

                dkda= dkda0 + dkda1 + dkda2
            endif
        endif

    end subroutine dKdAcalc

!+++-----------------------------------------------------
!+ Computation of g*dI2/dA,p81,RM1_MESH
!+++-----------------------------------------------------
!    subroutine gdI2dAcalc(ic, jc, yxs, Axs, gdI2dA)

!        implicit none

!        integer, intent(in) :: ic, jc
!        real, intent(in) :: yxs, Axs
!        real, intent(out) :: gdI2dA
!        real :: bwd, sslp, ybf

!        bwd= bo(ic,jc)
!        sslp= traps(ic,jc)

!        if (sslp==0.0) then
        !* rectangular channel
!            gdI2dA=grav*(Axs/(bwd**2.0))*dbdx(ic)
!        else
!            ybf=(Tw(ic,jc) - bwd)/(2.0*sslp)
!            if (yxs<=ybf) then
            !* inbank flow in trapezoidal channel
!                gdI2dA=grav*(1.0- bwd*(bwd**2.0 + 4.0*sslp*Axs)**(-0.5))/(2.0*sslp)*dbdx(ic)
!            else
            !* overbank flow in compound ch. with trap. main and rect. floodplains, p84-2,RM1
!                gdI2dA= grav*(dbdxCCm(ic)*ybf +  dbdxCCf(ic)*(yxs-ybf))/TwCC(ic,jc)
!            endif
!        endif

!    end subroutine gdI2dAcalc

!+++-----------------------------------------------------
!+ Computation of g*dI2/dA,p81,RM1_MESH
!+++-----------------------------------------------------
    subroutine gdI2dAcalc_mod(ic, yxs, Axs, dsgmdx, dsgmdxm, dsgmdxf, gdI2dA)

        implicit none

        integer, intent(in) :: ic
        real, intent(in) :: yxs, Axs, dsgmdx, dsgmdxm, dsgmdxf
        real, intent(out) :: gdI2dA
        real :: bwd, sslp, ybf

        bwd= bo(ic)
        sslp= traps(ic)

        if (sslp==0.0) then
        !* rectangular channel
            gdI2dA=grav*(Axs/(bwd**2.0))*dsgmdx !dbdx(ic)
        else
            ybf=(Tw(ic) - bwd)/(2.0*sslp)
            if (yxs<=ybf) then
            !* inbank flow in trapezoidal channel
                gdI2dA=grav*(1.0- bwd*(bwd**2.0 + 4.0*sslp*Axs)**(-0.5))/(2.0*sslp)*dsgmdx  !dbdx(ic)
            else
            !* overbank flow in compound ch. with trap. main and rect. floodplains, p84-2,RM1
                gdI2dA= grav*(dsgmdxm*ybf +  dsgmdxf*(yxs-ybf))/TwCC(ic)
            endif
        endif

    end subroutine gdI2dAcalc_mod

!+++------------------------------------------------------------
!+ Computation of Froud number of compound channel,p93~93-3,RM1
!+++------------------------------------------------------------
    subroutine FrCCcalc(ic, yxs, Axs, Qxs, FrCC)
    !** applied not only for compound channel having trapezoidal main channel with
    !*  rectangular floodplains, but also main channel with regular shapes
        implicit none

        integer, intent(in) :: ic
        real, intent(in) ::  yxs, Axs, Qxs
        real, intent(out) :: FrCC
        real :: bwd, sslp, Twd, TCCwd
        real :: ybf, TwCCi, A0, A1, A2, P0, P1, P2, K0, K1, K2, tlA, tlK
        real :: beta, betap, dA0dh, dK0dh, dA1dh, dK1dh, dA2dh, dK2dh, dAdh, dKdh
        real :: KA0, KA1, KA2, trm1, trm2, Vxs, Bxs, Dxs

        bwd= bo(ic)
        sslp= traps(ic)
        Twd= Tw(ic)
        TCCwd= TwCC(ic)

        if (sslp==0.0) then
        !* inbank flow in rec. channel
            Vxs= Qxs/Axs
            Bxs= bwd  !*top water-surface width
            Dxs= Axs/Bxs
            FrCC= Vxs/(grav*Dxs)**0.5
        else
            ybf=(Twd - bwd)/(2.0*sslp)
            if (yxs<ybf) then
            !* inbank flow in trap. main channel
                Vxs= Qxs/Axs
                Bxs= bwd + 2.0*sslp*yxs !*top water-surface width
                Dxs= Axs/Bxs
                FrCC= Vxs/(grav*Dxs)**0.5
            else
            !* overbank flow
                !** 1) beta, p93 ~ 93-1,RM1
                !* A0
                A0 = (bwd + sslp*ybf)*ybf + (yxs - ybf)*Twd
                !* K0
                P0 = bwd + 2.0*ybf*(1.0 + sslp**2.0)**0.5
                K0 = sk(ic)*A0**(5.0/3.0)/P0**(2.0/3.0)
                !* A1
                TwCCi=(TCCwd - Twd)/2.0
                A1=(yxs - ybf)*TwCCi
                !* K1, conveyance in the left floodplain (assuming symmetric)
                P1=TwCCi +  (yxs - ybf)
                K1 = skCC1(ic)*A1**(5.0/3.0)/P1**(2.0/3.0)
                !* A2
                TwCCi=(TCCwd - Twd)/2.0
                A2=(yxs - ybf)*TwCCi
                !* K2, conveyance in the right floodplain (assuming symmetric)
                P2=TwCCi +  (yxs - ybf)
                K2 = skCC2(ic)*A2**(5.0/3.0)/P2**(2.0/3.0)

                tlA= A0+A1+A2
                tlK= K0+K1+K2
                beta=(K0**2.0/A0 + K1**2.0/A1 + K2**2.0/A2)*tlA/(tlK**2.0)
                !** 2) beta prime, p93-1~-2,RM1
                dA0dh= Twd
                dK0dh= sk(ic)*(5.0/3.0)*Twd*(A0/P0)**(2.0/3.0)

                TwCCi= (TCCwd - Twd)/2.0
                dA1dh= TwCCi
                dK1dh = skCC1(ic)*((5.0/3.0)*TwCCi*(A1/P1)**(2.0/3.0) - (2.0/3.0)*(A1/P1)**(5.0/3.0))

                TwCCi= (TCCwd - Twd)/2.0
                dA2dh= TwCCi
                dK2dh= skCC2(ic)*((5.0/3.0)*TwCCi*(A2/P2)**(2.0/3.0) - (2.0/3.0)*(A2/P2)**(5.0/3.0))

                dAdh= dA0dh + dA1dh + dA2dh
                dKdh= dK0dh + dK1dh + dK2dh

                KA0= K0/A0
                KA1= K1/A1
                KA2= K2/A2

                trm1= (KA0*(2.0*dK0dh - KA0*dA0dh) + KA1*(2.0*dK1dh - KA1*dA1dh) + KA2*(2.0*dK2dh - KA2*dA2dh))*tlA/(tlK**2.0)
                trm2= (K0**2.0/A0 + K1**2.0/A1 + K2**2.0/A2)*(dAdh/(tlK**2.0) - 2.0*dKdh*tlA/(tlK**3.0))
                betap= trm1 + trm2
                Vxs= Qxs/Axs
                !** 3) Fc, p93,RM1
                FrCC= beta*Vxs/(grav*Axs/TCCwd + (beta*(beta-1) + Axs*betap/TCCwd )*Vxs**2.0)**0.5
            endif
        endif

    end subroutine FrCCcalc
 end module
