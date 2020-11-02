!Written  by Md Nazmul Azim Beg*, Ehab A Meselhe**
!*:     Postdoctoral Fellow, Tulane River and Coastal Center,
!       Department of River-Coastal Science and Engineering,
!       Tulane University
!**:    Professor, Department of River-Coastal Science and Engineering,
!       Tulane University

!Modified by Dong Ha Kim, NOAA's Office of Water Prediction, National Water Center
module diff

    use var
    use attrTable
    use subtools
    use nrtype

    implicit none

contains
!*--------------------------------------------------------------------------------
!*          Compute eei, ffi, exi, and fxi from upstream to downstream
!
!       input: dtini, dx, celty, diffusivity, q, qpx, qlat
!
! %%%%% qlat takes a unit of m^2/sec
!*--------------------------------------------------------------------------------
    subroutine ef_calc
        implicit none

        integer(KIND=i4b) :: i
        real(KIND=dp) :: cour !,  theta
        real(KIND=dp) :: a1, a2, a3, a4, b1, b2, b3, b4, dd1, dd2, dd3, dd4, h1, h2, h3, h4
        real(KIND=dp) :: qy, qxy, qxxy, qxxxy, ppi, qqi, rri, ssi, sxi, cour2
        real(KIND=dp), dimension(:), allocatable :: eei, ffi, exi, fxi

        !open(unit=1, file="./output/ef_calc_fortran.txt", status='unknown')

        allocate(eei(ncomp), ffi(ncomp), exi(ncomp), fxi(ncomp))
        eei(1) = 1.0
        ffi(1) = 0.
        exi(1) = 0.
        fxi(1) = 0.
        do i = 2,ncomp
            !* Calculation a1...a4, up to h4
            cour = dtini / dx(i-1)
            cour2= abs( celty(i) ) * cour

            a1 = 3.0 * cour2 ** 2.0 - 2.0 * cour2 ** 3.0
            a2 = 1 - a1
            a3 = ( cour2 ** 2.0 - cour2 ** 3.0 ) * dx(i-1)
            a4 = ( -1.0 * cour2 + 2.0 * cour2 ** 2.0 - cour2 ** 3.0 ) * dx(i-1)

            b1 = ( 6.0 * cour2 - 6.0 * cour2 ** 2.0 ) / ( -1.0 * dx(i-1) )
            b2 = - b1
            b3 = ( 2.0 * cour2 - 3.0 * cour2 ** 2.0 ) * ( -1.0 )
            b4 = ( -1.0 + 4.0 * cour2 - 3.0 * cour2 ** 2.0 ) * ( -1.0 )

            dd1 = ( 6.0 - 12.0 * cour2 ) / ( dx(i-1) ** 2.0 )
            dd2 = - dd1
            dd3 = ( 2.0 - 6.0 * cour2 ) / dx(i-1)
            dd4 = ( 4.0 - 6.0 * cour2 ) / dx(i-1)

            h1 = 12.0 / ( dx(i-1) ** 3.0 )
            h2 = - h1
            h3 = 6.0 / ( dx(i-1) ** 2.0 )
            h4 = h3

            qy   = a1 * q(i-1) + a2 * q(i) + a3 * qpx(i-1) + a4 * qpx(i)
            qxy  = b1 * q(i-1) + b2 * q(i) + b3 * qpx(i-1) + b4 * qpx(i)
            qxxy = dd1* q(i-1) + dd2* q(i) + dd3* qpx(i-1) + dd4* qpx(i)
            qxxxy= h1 * q(i-1) + h2 * q(i) + h3 * qpx(i-1) + h4 * qpx(i)

            ppi = - theta * diffty(i) * dtini / ( dx(i-1) ** 2.0 )
            qqi = 1.0 - 2.0 * ppi
            rri = ppi

            ssi = qy  + dtini * diffty(i) * ( 1.0 - theta ) * qxxy + dtini * celty(i) * qlatj(i-1)
            sxi = qxy + dtini * diffty(i) * ( 1.0 - theta ) * qxxxy+ dtini * celty(i) * qlatj(i-1)/ dx(i-1)

            eei(i) = -1.0 * rri / ( ppi * eei(i-1) + qqi )                     !! copied from split operator method
            ffi(i) = ( ssi - ppi * ffi(i-1) ) / ( ppi * eei(i-1) + qqi )       !! copied from split operator method
            exi(i) = -1.0 * rri / ( ppi * exi(i-1) + qqi )
            fxi(i) = ( sxi - ppi * fxi(i-1) ) / ( ppi * exi(i-1) + qqi )
        enddo

        !*----------------------------------------------------------------------------------
        !* compute q and qpx at ts+1 from downstream to upstream nodes for a single j reach
        !*----------------------------------------------------------------------------------
        q(ncomp)=q(ncomp-1) + qlatf*dx(ncomp-1) !* qlatf at ts+1 in unit of m2/sec
        qpx(ncomp)=0.0
        do i=ncomp-1, 1, -1
            q(i) = eei(i) * q(i+1) + ffi(i)    !* Q(i) at time t+dtini/60 [min] (=n+1)
            qpx(i)= exi(i) *qpx(i+1) + fxi(i)  !* Qx(i) at time t+dtini/60 [min] (=n+1)
        enddo

        deallocate(eei, ffi, exi, fxi)

    end subroutine ef_calc
!*--------------------------------------------------------------------------------
!*      Compute water elevation from downstream to upstream
!
!       INPUT: elv(ncomp), q, sk, dx
!       OUTPUT: elv(1 ~ ncomp-1), celty, diffty
!*--------------------------------------------------------------------------------
    subroutine elv_calc
        implicit none
        integer(KIND=i4b) :: i, icol, i1
        real(KIND=dp) ::  q_sk_multi, sfi, xt , width
        real(KIND=dp) :: depth_tz, area_tz, hydr_tz, conv_tz, hbf_tz, ybf_tz, topwd_tz, emann_tz, sum1, sum2, celav, diffav
        real(KIND=dp), allocatable ::  co(:)
        real(KIND=dp) :: So, hnm0, dsc, hnm, dmyq
        real(KIND=dp), dimension(:), allocatable ::  harr_m, qarr_m, harr_f, qarr_f
        real(KIND=dp), dimension(:), allocatable ::  elvarr, qarr

        allocate(co(ncomp))
        allocate(harr_m(nhincr_m), qarr_m(nhincr_m), harr_f(nhincr_f), qarr_f(nhincr_f))
        allocate(elvarr(nel),  qarr(nel))

        !open(unit=11, file="./output/elv_calc1_fortran.txt", status="unknown")
        !open(unit=12, file="./output/elv_calc2_fortran.txt", status="unknown")

        !y_opt=1 !* 1 for normal depth(kinematic); 2 for dept of diffusive wave.
        sum1=0.0
        sum2=0.0
        q_sk_multi = 1.0

        do i=ncomp,2,-1

            xt= elv(i)

            z0= z_ar(i)
            bo0= bo_ar(i)
            traps0= traps_ar(i)
            tw0= tw_ar(i)
            twcc0= twcc_ar(i)
            mann= mann_ar(i)
            manncc= manncc_ar(i)
            hbf_tz= (tw0 - bo0)/(2.0*traps0) !* bankfull depth

            if (tzeq_flag==0) then
            !* use x-sec attribute lookup tables
                icol= 2 !<-5
                co(i)= q_sk_multi*intp_elev(icol, i, xt)
                icol= 3 !<-6
                width= intp_elev(icol, i, xt)

            elseif (tzeq_flag==1) then
            !* use trapezoidal x-sec equations
                depth_tz= xt - z0
                call areacalc(depth_tz, area_tz)
                call hydRcalc(depth_tz, area_tz, hydr_tz)
                call Kcalc(depth_tz, area_tz, hydr_tz, conv_tz)
                if (depth_tz>hbf_tz) then
                    topwd_tz= twcc0
                else
                    topwd_tz= bo0 + 2.0*traps0*depth_tz
                end if
                co(i)= q_sk_multi*conv_tz
                width=topwd_tz
            endif

            sfi = (q(i) / co(i)) ** 2.0
            ybf_tz= hbf_tz + z_ar(i) !* bankfull depth/elevation for trapezoidal main channel
            if (xt.le.ybf_tz) then
            !* use Manning's N for main channel
                celty(i)=5.0 / 3.0 * sfi ** 0.3 * abs(q(i)) ** 0.4 / width ** 0.4 / ( mann*(1/q_sk_multi)) ** 0.6
            else
            !* when overbankflow, use equivalent manning's N
                dsc= q(i)
                depth_tz= xt - z_ar(i)
                emann_tz= emann_tmrf(depth_tz)
                celty(i)=5.0 / 3.0 * sfi ** 0.3 * abs(q(i)) ** 0.4 / width ** 0.4 / ( emann_tz*(1/q_sk_multi)) ** 0.6
            end if
            sum1= sum1 + celty(i)

            diffty(i) = abs(q(i)) / 2.0 / width / sfi
            sum2= sum2+ diffty(i)

            !* Y(i) at time t+dtini/60 [min] (=n+1)
            if (y_opt==1) then
            !* run normal depth computation at node i-1
                dsc= q(i-1)
                !* 2) normal depth by lookup tables
                if (tzeq_flag==0) then
                !* lookup table created from x-section approximation.
                    do i1=1, nel
                        elvarr(i1)= xsec_attr_rch(1, i1, i-1) !* elevation  i1=1,...,nel
                        qarr(i1)= xsec_attr_rch(4, i1, i-1)  !* uniform flow
                    enddo
                    elv(i-1)= intp_y(nel, qarr, elvarr, dsc)

                elseif (tzeq_flag==1) then
                !* lookup table created from x-section equations for trap.main & rect.floodplains
                    if (dsc<= ufqlt_m(nhincr_m, i-1)) then
                    !* inbank flow
                        do i1=1,nhincr_m
                            qarr_m(i1)=  ufqlt_m(i1, i-1)
                            harr_m(i1)= ufhlt_m(i1, i-1)
                        enddo
                        hnm= intp_y(nhincr_m, qarr_m, harr_m, dsc)
                    else
                    !* overbank flow
                        do i1=1,nhincr_f
                            qarr_f(i1)= ufqlt_f(i1, i-1)
                            harr_f(i1)= ufhlt_f(i1, i-1)
                        enddo
                        hnm= intp_y(nhincr_f, qarr_f, harr_f, dsc)
                    endif
                    elv(i-1)=hnm + z_ar(i-1)

                endif
            elseif (y_opt==2) then
            !* run diffusive depth computation
                elv(i-1) = elv(i) + sign ( sfi, q(i) ) * dx(i-1)
            endif
        enddo   !*do i=ncomp,2,-1

        !* averaged celty thorough all segments to prevent a wavefront of a cell with too high speed which
        !* produces wave breaking and instability
        celav= sum1/real(ncomp-1,KIND(sum1))
        diffav= sum2/real(ncomp-1,KIND(sum2))
        do i=1,ncomp
            celty(i) = celav
            diffty(i)= diffav
            if (diffty(i)>10.0) diffty(i)=10.0
        end do

        deallocate(co)
        deallocate(harr_m, qarr_m, harr_f, qarr_f)

    end subroutine elv_calc
!*----------------------------------------------------
!           Interpolate channel property by elev
!
!*----------------------------------------------------
    real function intp_elev(icol, inode, x)
        implicit none
        integer(KIND=i4b), intent(in) :: icol, inode
        real(KIND=dp), intent(in) :: x
        real(KIND=dp), allocatable :: dmyv(:)
        integer(KIND=i4b) :: i, i1, irow, nrow
        real(KIND=dp) :: x1, y1, x2, y2, y

        nrow= nel
        allocate(dmyv(nrow))
        do i1=1, nrow
            dmyv(i1)= xsec_attr_rch(1,i1,inode) !* elevation vector
        end do
        irow= locate(dmyv, x)
        if (irow.eq.0) irow= 1
        if (irow.eq.nrow) irow= nrow-1
        x1= xsec_attr_rch(1, irow, inode)
        y1= xsec_attr_rch(icol, irow, inode)
        x2= xsec_attr_rch(1, irow+1, inode)
        y2= xsec_attr_rch(icol, irow+1, inode)
        y= LInterpol(x1,y1,x2,y2,x)
        intp_elev= y
        deallocate(dmyv)

    end function intp_elev

end module diff
