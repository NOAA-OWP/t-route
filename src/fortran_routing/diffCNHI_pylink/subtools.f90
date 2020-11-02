!Written by Dong Ha Kim, NOAA's Office of Water Prediction, National Water Center

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
    function locate(xx,x)
        implicit none
        real(KIND=dp), dimension(:), intent(in) :: xx
        real(KIND=dp), intent(in) :: x
        integer(KIND=i4b) :: locate
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
        integer(KIND=i4b) :: i, irow
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

        bwd= bo0 !*bottom width
        sslp= traps0 !*trapezoidal
        tw= tw0
        twcc=twcc0

        if (sslp==0.0) then
        !* rectangular channel
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

        bxs=bo0
        ssxs=traps0
        twxs= tw0
        twccxs= twcc0

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

        bwd= bo0
        sslp= traps0
        Twd=tw0 !* top width of main channel
        TCCwd=twcc0 !* top width of compound channel
        hbf= (Twd - bwd)/(2.0*sslp)!* bankfull hxs

        if ((sslp==0.0).or.(hxs<=hbf)) then
        !* inbank flow in rectangular or trapezoidal channel
            cnvey=(1.0/mann)*Axs*(hydR**(2.0/3.0))
        else
        !* overbank flow in compound channel having trapezoidal main channel with rectangular floodplains, p84-2~p84-2-1, RM1
            !* conveyance in the main channel, K0
            subA = (bwd + sslp*hbf)*hbf + (hxs - hbf)*Twd
            subP = bwd + 2.0*hbf*(1.0+ sslp**2.0)**0.5
            K0 = (1.0/mann)*subA**(5.0/3.0)/subP**(2.0/3.0)
            !* conveyance in the left floodplain (assuming symmetric), K1
            TwCCi=(TCCwd - Twd)/2.0
            subA=(hxs - hbf)*TwCCi
            subP=TwCCi +  (hxs - hbf)
            K1 = (1.0/manncc)*subA**(5.0/3.0)/subP**(2.0/3.0)
            !* conveyance in the left floodplain (assuming symmetric), K1
            TwCCi=(TCCwd - Twd)/2.0
            subA=(hxs - hbf)*TwCCi
            subP=TwCCi +  (hxs - hbf)
            K2 = (1.0/manncc)*subA**(5.0/3.0)/subP**(2.0/3.0)

            cnvey=K0+K1+K2
        endif

    end subroutine Kcalc
!*----------------------------------------------------------------------
!*          Uniform flow discharge for trapz. main and rect. floodplain

!*  p.102, RM4
!*----------------------------------------------------------------------
    subroutine ufQ_tmrf(hxs, ufQ)
        implicit none
        real(KIND=dp), intent(in) :: hxs
        real(KIND=dp), intent(out) :: ufQ
        real(KIND=dp) :: Axs, hydR, WP, Twd, TCCwd, So
        real(KIND=dp) :: bwd, sslp, hbf, ufQ1, ufQ2, ufQ3

        bwd= bo0
        sslp= traps0
        Twd= tw0  !* top width of main channel
        TCCwd= twcc0 !* top width of compound channel
        So= So0
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
            ufQ2= (1.0/mann)*Axs*(hydR**(2.0/3.0))*(So**0.5)
        else
        !*overbank flow on rect. floodplains
            !* subsection 1 and 3, p102,RM4
            Axs= (hxs-hbf)*(TCCwd -Twd)/2.0
            WP= (TCCwd - Twd)/2.0 + (hxs-hbf)
            hydR= Axs/WP
            ufQ1= (1.0/manncc)*Axs*(hydR**(2.0/3.0))*(So**0.5)
            ufQ3= ufQ1
            !* subsection 2, p102,RM4
            Axs= (bwd + sslp*hbf)*hbf + (hxs-hbf)*Twd
            WP= bwd + 2.0*hbf*((1.0 + sslp**2.0)**0.5)
            hydR= Axs/WP
            ufQ2= (1.0/mann)*Axs*(hydR**(2.0/3.0))*(So**0.5)
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
        real(KIND=dp), intent(in) :: hxs
        real(KIND=dp) :: bwd, sslp, Twd, TCCwd
        real(KIND=dp) :: hbf, TwCCi, P0,P1,P2, tlP
        real(KIND=dp) ::  nom

        bwd= bo0
        sslp= traps0
        Twd= tw0   !* top width of main channel
        TCCwd= twcc0 !* top width of compound channel
        !* bankfull depth
        hbf=(Twd - bwd)/(2.0*sslp)
        !*P0,P1,P2 for overbank flow p84-2-1,RM1
        P0 = bwd + 2.0*hbf*(1.0 + sslp**2.0)**0.5
        TwCCi= (TCCwd-Twd)/2.0
        P1=TwCCi +  (hxs - hbf)
        P2=TwCCi +  (hxs - hbf)

        tlP= P0+P1+P2
        !* equivalent Manning's N for compound channel, Eq(4-35),p91,Chaudhry
        nom= P0*mann**(3.0/2.0) + P1*manncc**(3.0/2.0) + P2*manncc**(3.0/2.0)
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

    end function input_crtor1

end module subtools
