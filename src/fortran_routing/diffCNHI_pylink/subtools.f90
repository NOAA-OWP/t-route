MODULE subtools
    use var
    use nrtype

    implicit none

CONTAINS
!*-----------------------------------------------------------------------------
!               Locate function in f90, p.1045,NR f90
!
!   klo=max(min(locate(xa,x),n-1),1) In the Fortran 77 version of splint,
!   there is in-line code to find the location in the table by bisection. Here
!   we prefer an explicit call to locate, which performs the bisection. On
!   some massively multiprocessor (MMP) machines, one might substitute a different,
!   more parallel algorithm (see next note).
!*-----------------------------------------------------------------------------
    integer(KIND=i4b) function locate(xx,x)
        implicit none
        real(KIND=dp), dimension(:), intent(in) :: xx
        real(KIND=dp), intent(in) :: x
        !* Given an array xx(1:N), and given a value x, returns a value j such that x is between
        !* xx(j) and xx(j + 1). xx must be monotonic, either increasing or decreasing.
        !* j = 0 or j = N is returned to indicate that x is out of range.
        integer(KIND=i4b) :: n,jl,jm,ju
        logical :: ascnd

        n=size(xx)
        ascnd = (xx(n) >= xx(1))  !* True if ascending order of table, false otherwise.
        jl=0    !* Initialize lower
        ju=n+1  !* and upper limits.
        do
            if (ju-jl <= 1) exit    !* Repeat until this condition is satisfied.
            jm=(ju+jl)/2            !* Compute a midpoint,
            if (ascnd .eqv. (x >= xx(jm))) then
                jl=jm               !* and replace either the lower limit
            else
                ju=jm               !* or the upper limit, as appropriate.
            end if
        end do

        if (x == xx(1)) then        !* Then set the output, being careful with the endpoints.
            locate=1
        else if (x == xx(n)) then
            locate=n-1
        else
            locate=jl
        end if

    end function locate
!*--------------------------------------------------
!*      Same as function locate but subroutine
!
!*--------------------------------------------------
!    subroutine locate2(xx,x)
!        use nrtype
!        implicit none
!        real(KIND=dp)(SP), dimension(:), intent(in) :: xx
!        real(KIND=dp)(SP), intent(in) :: x
!        integer(KIND=i4b)(I4B) :: locate
!        !* Given an array xx(1:N), and given a value x, returns a value j such that x is between
!        !* xx(j) and xx(j + 1). xx must be monotonic, either increasing or decreasing.
!        !* j = 0 or j = N is returned to indicate that x is out of range.
!        integer(KIND=i4b)(I4B) :: n,jl,jm,ju
!        logical :: ascnd
!
!        n=size(xx)
!        ascnd = (xx(n) >= xx(1))  !* True if ascending order of table, false otherwise.
!        jl=0    !* Initialize lower
!        ju=n+1  !* and upper limits.
!        do
!            if (ju-jl <= 1) exit    !* Repeat until this condition is satisfied.
!            jm=(ju+jl)/2            !* Compute a midpoint,
!            if (ascnd .eqv. (x >= xx(jm))) then
!                jl=jm               !* and replace either the lower limit
!            else
!                ju=jm               !* or the upper limit, as appropriate.
!            end if
!        end do
!
!        if (x == xx(1)) then        !* Then set the output, being careful with the endpoints.
!            locate=1
!        else if (x == xx(n)) then
!            locate=n-1
!        else
!            locate=jl
!        end if
!        irow= locate  !* irow declared in var_module
!    end subroutine locate2
!*--------------------------------------------------
!*                 Linear Interpolation
!
!*--------------------------------------------------
    real(KIND=dp) function LInterpol(x1,y1,x2,y2,x)
        implicit none
        real(KIND=dp), intent(in) :: x1, y1, x2, y2, x

        !* interpolate y for the given x
        LInterpol= (y2-y1)/(x2-x1)*(x-x1)+y1

    end function LInterpol
!*--------------------------------------------
!           Interpolate any value
!
!*--------------------------------------------
    real(KIND=dp) function intp_y(nrow, xarr, yarr, x)
        implicit none
        integer(KIND=i4b), intent(in) :: nrow
        real(KIND=dp), dimension(nrow), intent(in) :: xarr, yarr
        real(KIND=dp), intent(in) :: x
        integer(KIND=i4b) :: irow
        real(KIND=dp) :: x1, y1, x2, y2, y

        irow= locate(xarr, x)
        if (irow.eq.0) irow= 1
        if (irow.eq.nrow) irow= nrow-1
        x1= xarr(irow); y1= yarr(irow)
        x2= xarr(irow+1); y2= yarr(irow+1)
        y= LInterpol(x1,y1,x2,y2,x)
        intp_y = y

    end function intp_y
!+++---------------------------------------------------------------------------
!+ Computation of area of various channel x-sections with depth as an argument
!+ and without pre-specified chshp
!+++---------------------------------------------------------------------------
    subroutine areacalc(hxs, Axs)

        implicit none

        real(KIND=dp), intent(in) :: hxs
        real(KIND=dp), intent(out) :: Axs
        real(KIND=dp) :: bwd, sslp, hbf, tw, twcc

        bwd= bo_g !*bottom width
        sslp= traps_g !*trapezoidal
        tw= tw_g
        twcc=twcc_g

        if (sslp==0.0) then
        !* rectangular channel_g
            Axs=bwd*hxs
        else
        !* determine whether inbank or overbank flow
            !* bankfull depth
            hbf= (tw - bwd)/(2.0*sslp)
            if (hxs.le.hbf) then
                !* trapezoidal channel inbank flow
                Axs=(bwd + sslp*hxs)*hxs
            else
                !*overbank flow on rect. floodplains
                Axs=(bwd + sslp*hbf)*hbf + twcc*(hxs-hbf)
            end if
        endif

    end subroutine areacalc
!+++-----------------------------------------------------
!+ Computation of hydraulic radius R (=A/P)
!+++-----------------------------------------------------
    subroutine hydRcalc(hxs, Axs, hydR)
        implicit none
        real(KIND=dp), intent(in) ::  hxs, Axs
        real(KIND=dp), intent(out) :: hydR
        real(KIND=dp) :: bxs, ssxs, twxs, twccxs
        real(KIND=dp) :: hbf

        bxs=bo_g
        ssxs=traps_g
        twxs= tw_g
        twccxs= twcc_g

        if (ssxs==0.0) then
        !* rectangular channel
            hydR= Axs/(bxs + 2.0*hxs)
        else
            hbf= (twxs - bxs)/(2.0*ssxs)
            if (hxs<=hbf) then
            !* inbank flow in trapezoidal channel
                hydR=Axs/(bxs + 2.0*hxs*((1.0 + ssxs**2.0)**0.5))
            else
            !* compound channel having trapezoidal main and rectangular floodplains, p80-,RM1
                hydR= Axs/(bxs + 2.0*hbf*(1.0 + ssxs**2.0)**0.5 + twccxs - twxs + 2.0*(hxs-hbf))
            endif
        endif

    end subroutine hydRcalc
!+++-----------------------------------------------------
!+ Computation of conveyance K (=1/N * A * R^(2/3))
!+++-----------------------------------------------------
    subroutine Kcalc(hxs, Axs, hydR, cnvey)

        implicit none
        real(KIND=dp), intent(in) :: hxs, Axs, hydR
        real(KIND=dp), intent(out) :: cnvey
        real(KIND=dp) :: bwd, sslp, Twd, TCCwd
        real(KIND=dp) :: subA, subP, hbf, TwCCi, K0, K1, K2

        bwd= bo_g
        sslp= traps_g
        Twd=tw_g !* top width of main channel
        TCCwd=twcc_g !* top width of compound channel
        hbf= (Twd - bwd)/(2.0*sslp)!* bankfull hxs

        if ((sslp==0.0).or.(hxs<=hbf)) then
        !* inbank flow in rectangular or trapezoidal channel
            cnvey=(1.0/mann_g)*Axs*(hydR**(2.0/3.0))
        else
        !* overbank flow in compound channel having trapezoidal main channel with rectangular floodplains, p84-2~p84-2-1, RM1
            !* conveyance in the main channel, K0
            subA = (bwd + sslp*hbf)*hbf + (hxs - hbf)*Twd
            subP = bwd + 2.0*hbf*(1.0+ sslp**2.0)**0.5
            K0 = (1.0/mann_g)*subA**(5.0/3.0)/subP**(2.0/3.0)
            !* conveyance in the left floodplain (assuming symmetric), K1
            TwCCi=(TCCwd - Twd)/2.0
            subA=(hxs - hbf)*TwCCi
            subP=TwCCi +  (hxs - hbf)
            K1 = (1.0/manncc_g)*subA**(5.0/3.0)/subP**(2.0/3.0)
            !* conveyance in the left floodplain (assuming symmetric), K1
            TwCCi=(TCCwd - Twd)/2.0
            subA=(hxs - hbf)*TwCCi
            subP=TwCCi +  (hxs - hbf)
            K2 = (1.0/manncc_g)*subA**(5.0/3.0)/subP**(2.0/3.0)

            cnvey=K0+K1+K2
        endif

    end subroutine Kcalc


!+++-------------------------------------------------------------------------------
!+ computation of normal depth in regular/trapezoidal/compound channel using
!+ Newton-Raphson/Approximated/Bisection methods.
!+ Refer to Appendix C-2, Chaudhary and p71,RM1_MESH
!+++-------------------------------------------------------------------------------
    !subroutine nmaldepcalc(ic, jc, So, ynm0, dsc, ynm)
    subroutine nmaldepcalc(So, ynm0, dsc, ynm)
        implicit none

        !integer(KIND=i4b), intent(in) :: ic, jc
        real(KIND=dp), intent(in) :: So, ynm0, dsc
        real(KIND=dp), intent(out) :: ynm

        integer(KIND=i4b) :: recnm_app, ib, ib_mx
        real(KIND=dp) :: sslp, bwd !, manN
        real(KIND=dp) :: c0, c1
        real(KIND=dp) :: errNR, tolNR
        real(KIND=dp) :: lambda, epsr, zetn
        real(KIND=dp) :: yb0, yb1, yb2, fyn, yinc, ybp, mpy, cvrb, errb


        tolNR=0.01  !*tolerance
        c0=1.0    !*1 for SI units and 1.49 for customary English units
        lambda= c0
    !manN= mann !1/sk(ic,jc)
        !So= -(z(ic,jc) - z(ic-1,jc))/dx_ar_g(ic-1,jc) !*this subroutine is called only in predictor step ** changed 5/13/2020
        !So= chbtslp(ic-1,jc)
        sslp= traps_g !traps(ic,jc)
        bwd= bo_g !bo(ic,jc)
        c1=mann_g*dsc/(c0*(So**0.5))
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
                epsr=mann_g*dsc/(lambda*(bwd**(8.0/3.0))*(So**0.5))                         !*Eq.8
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
            fyn= fCCyn(So, ynm0, dsc) !fCCyn(ic, jc, So, ynm0, dsc)
            if (fyn==0.0) then
                ynm=ynm0
                return
            elseif (fyn>0.0) then
                yb2=ynm0
                !* search for yb1
                yinc=0.1*ynm0
                ybp= ynm0 - yinc
                fyn= fCCyn(So, ybp, dsc) !fCCyn(ic, jc, So, ybp, dsc)
                do while ((fyn>0.0).and.(ybp>yinc))
                    ybp=ybp - yinc
                    fyn= fCCyn(So, ybp, dsc) !fCCyn(ic, jc, So, ybp, dsc)
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
                fyn= fCCyn(So, ybp, dsc)  !fCCyn(ic, jc, So, ybp, dsc)
                mpy=10.0  !*multiplier for y for deciding possibly largest value of y
                do while ((fyn<0.0).and.(ybp<mpy*ynm0))
                    ybp=ybp+yinc
                    fyn= fCCyn(So, ybp, dsc)  ! fCCyn(ic, jc, So, ybp, dsc)
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
            ib_mx=40
            do ib=1,ib_mx
                yb0= (yb1 + yb2)/2.0
                fyn= fCCyn(So, yb0, dsc)  !fCCyn(ic, jc, So, yb0, dsc)
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
                elseif (ib==ib_mx) then
                !* when iter reaches its max,
                    ynm= (yb1 + yb2)/2.0
                    write(*,"(A50,2F10.4)") "norm.depth Bisection iter reaches max", ynm, fyn
                    !pause
                end if
            enddo

        endif

    end subroutine nmaldepcalc

!+++----------------------------------------------
!+ F(yn)=A**(5/3)/P**(2/3) - nQ/So**(1/2), p95,RM1
!+++----------------------------------------------
    !real(KIND=dp) function fCCyn(ic, jc, So, hxs, dsc)
    real(KIND=dp) function fCCyn(So, hxs, dsc)
        implicit none
        !integer(KIND=i4b), intent(in) :: ic, jc
        real(KIND=dp), intent(in) :: So, hxs, dsc
        real(KIND=dp) :: bwd, sslp, Twd, TCCwd
        real(KIND=dp) :: hbf, TwCCi, A0,A1,A2, P0,P1,P2, tlA, tlP
        real(KIND=dp) :: emanN, nom  !manN,


        bwd= bo_g ! bo(ic,jc)
        sslp= traps_g !traps(ic,jc)
        Twd= tw_g  !Tw(ic,jc) !* top width of main channel
        TCCwd= twcc_g !* top width of compound channel

        hbf=(Twd - bwd)/(2.0*sslp)
        if (hxs<=hbf) then
        !* still inbank flow
            !call areacalc(ic, jc, hxs, A0)
            call areacalc(hxs, A0)

            P0 = bwd + 2.0*hxs*(1.0 + sslp**2.0)**0.5
            !manN= 1.0/sk(ic,jc)
            fCCyn= A0**(5.0/3.0)/P0**(2.0/3.0) - mann_g*dsc/(So**0.5)
        else
        !* overbank flow
        !** A0,A1,A2 and P0,P1,P2 for overbank flow p84-2-1,RM1
            !* A0 & P0
            A0 = (bwd + sslp*hbf)*hbf + (hxs - hbf)*Twd
            P0 = bwd + 2.0*hbf*(1.0 + sslp**2.0)**0.5

            TwCCi= (TCCwd-Twd)/2.0 !(TwCC(ic,jc)-Twd)/2.0
            !* A1 & P1
            A1=(hxs - hbf)*TwCCi
            P1=TwCCi +  (hxs - hbf)
            !* A2 & P2
            A2= (hxs - hbf)*TwCCi
            P2=TwCCi +  (hxs - hbf)
            !** F(yn)
            tlA= A0+A1+A2
            tlP= P0+P1+P2
            !* equivalent Manning's N for compound channel, Eq(4-35),p91,Chaudhry
            !nom= P0*(1.0/sk(ic,jc))**(3.0/2.0) + P1*(1.0/skCC1(ic,jc))**(3.0/2.0) + P2*(1.0/skCC2(ic,jc))**(3.0/2.0)
            nom= P0*mann_g**(3.0/2.0) + P1*manncc_g**(3.0/2.0) + P2*manncc_g**(3.0/2.0)
            emanN= (nom/tlP)**(2.0/3.0)
            fCCyn= tlA**(5.0/3.0)/tlP**(2.0/3.0) - emanN*dsc/(So**0.5)
        endif

    end function fCCyn

!*----------------------------------------------------------------------
!*          Uniform flow discharge for trapz. main and rect. floodplain

!*  p.102, RM4
!*----------------------------------------------------------------------
    subroutine ufQ_tmrf(hxs, ufQ)
        implicit none
        !integer(KIND=i4b), intent(in) :: ic, jc
        real(KIND=dp), intent(in) :: hxs
        real(KIND=dp), intent(out) :: ufQ
        real(KIND=dp) :: Axs, hydR, WP, Twd, TCCwd, So
        real(KIND=dp) :: bwd, sslp, hbf, ufQ1, ufQ2, ufQ3

        bwd= bo_g ! bo(ic,jc)
        sslp= traps_g !traps(ic,jc)
        Twd= tw_g  !Tw(ic,jc) !* top width of main channel
        TCCwd= twcc_g !* top width of compound channel
        So= So_g
        !mann= mna(i,j)
        !manncc= mncca(i,j)

        !* determine whether inbank or overbank flow
        !* bankfull depth
        hbf= (Twd - bwd)/(2.0*sslp)
        if (hxs.le.hbf) then
        !* trapezoidal channel inbank flow
            Axs=(bwd + sslp*hxs)*hxs
            WP= bwd + 2.0*hxs*((1.0+(sslp**2.0))**0.5)
            hydR=Axs/WP
            ufQ1=0.0
            ufQ3=0.0
            ufQ2= (1.0/mann_g)*Axs*(hydR**(2.0/3.0))*(So**0.5)
        else
        !*overbank flow on rect. floodplains
            !* subsection 1 and 3, p102,RM4
            Axs= (hxs-hbf)*(TCCwd -Twd)/2.0
            WP= (TCCwd - Twd)/2.0 + (hxs-hbf)
            hydR= Axs/WP
            ufQ1= (1.0/manncc_g)*Axs*(hydR**(2.0/3.0))*(So**0.5)
            ufQ3= ufQ1
            !* subsection 2, p102,RM4
            Axs= (bwd + sslp*hbf)*hbf + (hxs-hbf)*Twd
            WP= bwd + 2.0*hbf*((1.0 + sslp**2.0)**0.5)
            hydR= Axs/WP
            ufQ2= (1.0/mann_g)*Axs*(hydR**(2.0/3.0))*(So**0.5)
        end if

        ufQ= ufQ1+ufQ2+ufQ3

    end subroutine ufQ_tmrf

!*--------------------------------------------------------------------------
!*          Equivalent Manning's N for overbank flow in
!*                                      trapz.main and rect. floodplain

!*  Eq(4-35),p91,Chaudhry p95,RM1
!*--------------------------------------------------------------------------
    real(KIND=dp) function emann_tmrf(hxs)
        implicit none
        !integer(KIND=i4b), intent(in) :: ic, jc
        real(KIND=dp), intent(in) :: hxs
        real(KIND=dp) :: bwd, sslp, Twd, TCCwd
        real(KIND=dp) :: hbf, TwCCi, P0,P1,P2, tlP
        real(KIND=dp) ::  nom

        bwd= bo_g ! bo(ic,jc)
        sslp= traps_g !traps(ic,jc)
        Twd= tw_g  !Tw(ic,jc) !* top width of main channel
        TCCwd= twcc_g !* top width of compound channel
        !mann= mna(i,j)
        !manncc= mncca(i,j)

        !* bankfull depth
        hbf=(Twd - bwd)/(2.0*sslp)
        !*P0,P1,P2 for overbank flow p84-2-1,RM1
        P0 = bwd + 2.0*hbf*(1.0 + sslp**2.0)**0.5
        TwCCi= (TCCwd-Twd)/2.0
        P1=TwCCi +  (hxs - hbf)
        P2=TwCCi +  (hxs - hbf)

        tlP= P0+P1+P2
        !* equivalent Manning's N for compound channel, Eq(4-35),p91,Chaudhry
        nom= P0*mann_g**(3.0/2.0) + P1*manncc_g**(3.0/2.0) + P2*manncc_g**(3.0/2.0)
        emann_tmrf= (nom/tlP)**(2.0/3.0)

    end function emann_tmrf
!*--------------------------------------------------------------------------
!*        Artificial inflow or lateral or water depth data
!*--------------------------------------------------------------------------
    real(KIND=dp) function input_crtor1(mtp, mn, vr, yintc, n)
        implicit none
        integer(KIND=i4b), intent(in) :: n
        real(KIND=dp), intent(in) :: mtp, mn, vr, yintc
        real(KIND=dp) :: ifunc

        ifunc= mtp*exp(-0.5*((real(n,KIND(mn))-mn)/vr)**2.0)/(vr*(2.0*3.14)**0.5) + yintc
        input_crtor1= ifunc

        !qlat_g(n,i,j)= mtp*exp(-0.5*((real(n,KIND(mn))-mn)/vr)**2.0)/(vr*(2.0*3.14)**0.5) + 0.1
        !qlat_g(n,i,j)= qlat_g(n,i,j)/dx_ar_g(i,j)
        !j=1
        !qbd(n,1,j)= mtp*exp(-0.5*((real(n,KIND(mn))-mn)/vr)**2.0)/(vr*(2.0*3.14)**0.5)+2.0
        !j=2
        !qbd(n,1,j)= mtp*exp(-0.5*((real(n,KIND(mn))-mn)/vr)**2.0)/(vr*(2.0*3.14)**0.5)+4.0

        !ybd(n,ncomp,j) = mtp2*exp(-0.5*((real(n,KIND(mn2))-mn2)/vr2)**2.0)/(vr2*(2.0*3.14)**0.5)+1.0 + z(ncomp,j)

    end function input_crtor1

end module subtools
