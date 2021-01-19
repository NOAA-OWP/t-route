!            FINITE DIFFERENCE METHOD
!
!  A program for one dimensional flow in open channel
!
!program mesh
module diff


    use var
    !use arrays
    !use attrTable
    use subtools
    use nrtype

    implicit none

contains


!*--------------------------------------------------------------------------------
!*          Compute eei_g, ffi_g, exi_g, and fxi_g from upstream to downstream
!
!       input: dtini, dx_ar_g, celty_g, diffusivity, q, qpx_g, qlat_g
!
! %%%%% qlat_g takes a unit of m^2/sec
!*--------------------------------------------------------------------------------
    subroutine ef_calc(ncomp, j)
        implicit none

        integer(KIND=i4b),intent(in) :: ncomp, j

        integer(KIND=i4b) :: i
        real(KIND=dp) :: cour !,  theta
        real(KIND=dp) :: a1, a2, a3, a4, b1, b2, b3, b4, dd1, dd2, dd3, dd4, h1, h2, h3, h4
        real(KIND=dp) :: qy, qxy, qxxy, qxxxy, ppi, qqi, rri, ssi, sxi, cour2

        !open(unit=1, file="./output/ef_calc_fortran.txt", status='unknown')
        !open(unit=2, file="./output/ssi sxi.txt", status='unknown')
        !allocate(eei_g(ncomp), ffi_g(ncomp), exi_g(ncomp), fxi_g(ncomp))
                !write(1,"(2A20, A5, 8A20)") "pytime", "pyhsegid", "i",&
                !                    "ppi", "eei_g(i-1)", "exi_g(i-1)", "qqi",&
                !                    "eei_g(i)", "ffi_g(i)", "exi_g(i)", "fxi_g(i)"
                !write(2,"(2A20, A5, 12A20)") "pytime", "pyhsegid", "i",&
                !        "qy", "dtini", "diffty_g(i)", "theta", "qxxy", "celty_g(i)", "qlatj_g(i-1)", "ssi",&
                !        "qxy", "qxxxy", "dx_ar_g(i-1)", "sxi"
        !* initialize values
        !eei_g(1) = 1.0
        !ffi_g(1) = 0.
        !exi_g(1) = 0.
        !fxi_g(1) = 0.
        do i = 2,ncomp
            !* Calculation a1...a4, up to h4...
            cour = dtini_g / dx_ar_g(i-1,j)
            cour2= abs( celty_g(i,j) ) * cour

            a1 = 3.0 * cour2 ** 2.0 - 2.0 * cour2 ** 3.0
            a2 = 1 - a1
            a3 = ( cour2 ** 2.0 - cour2 ** 3.0 ) * dx_ar_g(i-1,j)
            a4 = ( -1.0 * cour2 + 2.0 * cour2 ** 2.0 - cour2 ** 3.0 ) * dx_ar_g(i-1,j)

            b1 = ( 6.0 * cour2 - 6.0 * cour2 ** 2.0 ) / ( -1.0 * dx_ar_g(i-1,j) )
            b2 = - b1
            b3 = ( 2.0 * cour2 - 3.0 * cour2 ** 2.0 ) * ( -1.0 )
            b4 = ( -1.0 + 4.0 * cour2 - 3.0 * cour2 ** 2.0 ) * ( -1.0 )

            dd1 = ( 6.0 - 12.0 * cour2 ) / ( dx_ar_g(i-1,j) ** 2.0 )
            dd2 = - dd1
            dd3 = ( 2.0 - 6.0 * cour2 ) / dx_ar_g(i-1,j)
            dd4 = ( 4.0 - 6.0 * cour2 ) / dx_ar_g(i-1,j)

            h1 = 12.0 / ( dx_ar_g(i-1,j) ** 3.0 )
            h2 = - h1
            h3 = 6.0 / ( dx_ar_g(i-1,j) ** 2.0 )
            h4 = h3

            qy   = a1 * q_g(i-1,j) + a2 * q_g(i,j) + a3 * qpx_g(i-1,j) + a4 * qpx_g(i,j)
            qxy  = b1 * q_g(i-1,j) + b2 * q_g(i,j) + b3 * qpx_g(i-1,j) + b4 * qpx_g(i,j)
            qxxy = dd1* q_g(i-1,j) + dd2* q_g(i,j) + dd3* qpx_g(i-1,j) + dd4* qpx_g(i,j)
            qxxxy= h1 * q_g(i-1,j) + h2 * q_g(i,j) + h3 * qpx_g(i-1,j) + h4 * qpx_g(i,j)

            ppi = - theta_g * diffty_g(i,j) * dtini_g / ( dx_ar_g(i-1,j) ** 2.0 )
            qqi = 1.0 - 2.0 * ppi
            rri = ppi

            ssi = qy  + dtini_g * diffty_g(i,j) * ( 1.0 - theta_g ) * qxxy + dtini_g * celty_g(i,j) * qlatj_g(i-1) !qlat_g(i,j)
            sxi = qxy + dtini_g * diffty_g(i,j) * ( 1.0 - theta_g ) * qxxxy+ dtini_g * celty_g(i,j) * qlatj_g(i-1)/ dx_ar_g(i-1,j)

            eei_g(i) = -1.0 * rri / ( ppi * eei_g(i-1) + qqi )                     !! copied from split operator method
            ffi_g(i) = ( ssi - ppi * ffi_g(i-1) ) / ( ppi * eei_g(i-1) + qqi )       !! copied from split operator method
            exi_g(i) = -1.0 * rri / ( ppi * exi_g(i-1) + qqi )
            fxi_g(i) = ( sxi - ppi * fxi_g(i-1) ) / ( ppi * exi_g(i-1) + qqi )
            !write(*, *) i, eei_g(i), ffi_g(i), exi_g(i), fxi_g(i)
        enddo

        !deallocate(eei_g, ffi_g, exi_g, fxi_g)

    end subroutine ef_calc

!*-------------------------------------------------------------------------------------
!*          Compute eei_m_g, ffi_m_g, exi_m_g, and fxi_m_g from upstream to downstream
!*          for modified segment lengths of a given reach j
!
!       input: dtini, dx_ar_m_g, celty_m_g, diffty_m_g, q_m_g, qpx_m_g, qlat_m_g
!
! %%%%% qlat_m_g takes a unit of m^2/sec
!*--------------------------------------------------------------------------------------
    subroutine ef_calc_m(ncomp_m)
        implicit none

        integer(KIND=i4b),intent(in) :: ncomp_m

        integer(KIND=i4b) :: i
        real(KIND=dp) :: cour !,  theta
        real(KIND=dp) :: a1, a2, a3, a4, b1, b2, b3, b4, dd1, dd2, dd3, dd4, h1, h2, h3, h4
        real(KIND=dp) :: qy, qxy, qxxy, qxxxy, ppi, qqi, rri, ssi, sxi, cour2

        !open(unit=1, file="./output/ef_calc_fortran.txt", status='unknown')
        !open(unit=2, file="./output/ssi sxi.txt", status='unknown')
        !allocate(eei_g(ncomp), ffi_g(ncomp), exi_g(ncomp), fxi_g(ncomp))
                !write(1,"(2A20, A5, 8A20)") "pytime", "pyhsegid", "i",&
                !                    "ppi", "eei_g(i-1)", "exi_g(i-1)", "qqi",&
                !                    "eei_g(i)", "ffi_g(i)", "exi_g(i)", "fxi_g(i)"
                !write(2,"(2A20, A5, 12A20)") "pytime", "pyhsegid", "i",&
                !        "qy", "dtini", "diffty_g(i)", "theta", "qxxy", "celty_g(i)", "qlatj_g(i-1)", "ssi",&
                !        "qxy", "qxxxy", "dx_ar_g(i-1)", "sxi"

        do i = 2, ncomp_m
            !* Calculation a1...a4, up to h4...
            cour = dtini_g / dx_ar_m_g(i-1)
            cour2= abs( celty_m_g(i) ) * cour

            a1 = 3.0 * cour2 ** 2.0 - 2.0 * cour2 ** 3.0
            a2 = 1 - a1
            a3 = ( cour2 ** 2.0 - cour2 ** 3.0 ) * dx_ar_m_g(i-1)
            a4 = ( -1.0 * cour2 + 2.0 * cour2 ** 2.0 - cour2 ** 3.0 ) * dx_ar_m_g(i-1)

            b1 = ( 6.0 * cour2 - 6.0 * cour2 ** 2.0 ) / ( -1.0 * dx_ar_m_g(i-1) )
            b2 = - b1
            b3 = ( 2.0 * cour2 - 3.0 * cour2 ** 2.0 ) * ( -1.0 )
            b4 = ( -1.0 + 4.0 * cour2 - 3.0 * cour2 ** 2.0 ) * ( -1.0 )

            dd1 = ( 6.0 - 12.0 * cour2 ) / ( dx_ar_m_g(i-1) ** 2.0 )
            dd2 = - dd1
            dd3 = ( 2.0 - 6.0 * cour2 ) / dx_ar_m_g(i-1)
            dd4 = ( 4.0 - 6.0 * cour2 ) / dx_ar_m_g(i-1)

            h1 = 12.0 / ( dx_ar_m_g(i-1) ** 3.0 )
            h2 = - h1
            h3 = 6.0 / ( dx_ar_m_g(i-1) ** 2.0 )
            h4 = h3

            qy   = a1 * q_m_g(i-1) + a2 * q_m_g(i) + a3 * qpx_m_g(i-1) + a4 * qpx_m_g(i)
            qxy  = b1 * q_m_g(i-1) + b2 * q_m_g(i) + b3 * qpx_m_g(i-1) + b4 * qpx_m_g(i)
            qxxy = dd1* q_m_g(i-1) + dd2* q_m_g(i) + dd3* qpx_m_g(i-1) + dd4* qpx_m_g(i)
            qxxxy= h1 * q_m_g(i-1) + h2 * q_m_g(i) + h3 * qpx_m_g(i-1) + h4 * qpx_m_g(i)

            ppi = - theta_g * diffty_m_g(i) * dtini_g / ( dx_ar_m_g(i-1) ** 2.0 )
            qqi = 1.0 - 2.0 * ppi
            rri = ppi

            ssi = qy  + dtini_g * diffty_m_g(i) * ( 1.0 - theta_g ) * qxxy + dtini_g * celty_m_g(i) * qlatj_m_g(i-1) !qlat_g(i,j)
            sxi = qxy + dtini_g * diffty_m_g(i) * ( 1.0 - theta_g ) * qxxxy+ dtini_g * celty_m_g(i) * qlatj_m_g(i-1)/ dx_ar_m_g(i-1)

            eei_m_g(i) = -1.0 * rri / ( ppi * eei_m_g(i-1) + qqi )                     !! copied from split operator method
            ffi_m_g(i) = ( ssi - ppi * ffi_m_g(i-1) ) / ( ppi * eei_m_g(i-1) + qqi )       !! copied from split operator method
            exi_m_g(i) = -1.0 * rri / ( ppi * exi_m_g(i-1) + qqi )
            fxi_m_g(i) = ( sxi - ppi * fxi_m_g(i-1) ) / ( ppi * exi_m_g(i-1) + qqi )
        enddo

    end subroutine ef_calc_m
!*--------------------------------------------------------------------------------
!*      Compute water elevation from downstream to upstream
!
!       INPUT: elv(ncomp), q, sk, dx_ar_g
!       OUTPUT: elv(1 ~ ncomp-1), celty_g, diffty_g
!*--------------------------------------------------------------------------------
    subroutine elv_calc(ts, ncomp, j)
        implicit none

        integer(KIND=i4b),intent(in) :: ts, ncomp, j

        integer(KIND=i4b) :: i, icol, i1
        real(KIND=dp) ::  q_sk_multi, sfi, xt , width
        real(KIND=dp) :: depth_tz, area_tz, hydr_tz, conv_tz, hbf_tz, ybf_tz, topwd_tz, emann_tz, sum1, sum2, celav, diffav
        real(KIND=dp), allocatable ::  co(:)
        real(KIND=dp) :: dsc, hnm

        real(KIND=dp), dimension(:), allocatable ::  harr_m, qarr_m, harr_f, qarr_f
        real(KIND=dp), dimension(:), allocatable ::  elvarr, qarr
        integer(KIND=i4b) :: j1, pncomp !for test

        allocate(co(ncomp))
        allocate(harr_m(nhincr_m_g), qarr_m(nhincr_m_g), harr_f(nhincr_f_g), qarr_f(nhincr_f_g))
        allocate(elvarr(nel_g),  qarr(nel_g))

!        open(unit=301, file="./output/elv_calc test.txt")


        !*-------------------------------------------------------------------------------
        !*    Partial diffusive elevation(=normal depth)
        !*    elv at i <- normal_depth_LT{q at i, ch.geo at i}, where LT is lookup table
        !*-------------------------------------------------------------------------------
!                write(301,"(A10,A2,A10,4A15)") "ts", "i", "j", "q_g(i,j)", "hnm", "z_ar_g(i,j)", "elv_g(i,j)"
        if (y_opt_g==1) then
            do i= ncomp-1, 1, -1 !* by p.3,RM5

                dsc= abs(q_g(i,j))

                !* By lookup tables
                if (tzeq_flag_g==0) then
                !* lookup table created from x-section approximation.
                    do i1=1, nel_g
                        elvarr(i1)= xsec_attr_rch_g(i, j, i1, 1)   !* elevation  i1=1,...,nel
                        qarr(i1)= xsec_attr_rch_g(i, j, i1, 4)     !* uniform flow
                    enddo
                    elv_g(i,j)= intp_y(nel_g, qarr, elvarr, dsc)

                elseif (tzeq_flag_g==1) then
                !* lookup table created from x-section equations for trap.main & rect.floodplains
                    if (dsc<= ufqlt_m_g(i, j, nhincr_m_g)) then
                    !* inbank flow
                        do i1=1,nhincr_m_g
                            qarr_m(i1)= ufqlt_m_g(i, j, i1)
                            harr_m(i1)= ufhlt_m_g(i, j, i1)
                        enddo
                        hnm= intp_y(nhincr_m_g, qarr_m, harr_m, dsc)
                    else
                        !* overbank flow
                        do i1=1,nhincr_f_g
                            qarr_f(i1)= ufqlt_f_g(i, j, i1)
                            harr_f(i1)= ufhlt_f_g(i, j, i1)
                        enddo
                        hnm= intp_y(nhincr_f_g, qarr_f, harr_f, dsc)
                    endif
                    elv_g(i,j)= hnm + z_ar_g(i,j)
                endif
            enddo !* do i= ncomp-1, 1, -1

!                    do i1=1,ncomp
!                        write(301,"(I10,I2,I10,4F15.3)") ts, i1, j, q_g(i1,j), hnm, z_ar_g(i1,j), elv_g(i1,j)
!                    enddo
        endif !* if (y_opt_g==1) then
!
!                   write(301,"(A10,A2,A10,10A15)") "ts", "i", "j", "q_g(i,j)", "xt", "depth_tz", "co(i)", "sfi", "width",&
!                                                "emann_tz", "celty_g(i,j)", "diffty_g(i,j)", "elv_g(i,j)"
        sum1=0.0
        sum2=0.0
        q_sk_multi = 1.0
        do i=ncomp-1,1,-1 !* by p.3,RM5

            if (y_opt_g==1) then
            !* normal depth computed above.
                xt= elv_g(i,j)
            elseif (y_opt_g==2) then
            !* for diffusive elevation method, elevation to compute Sf at i is approximated by elev. at i+1
                xt= elv_g(i+1,j)
            endif

            z_g= z_ar_g(i,j)
            bo_g= bo_ar_g(i,j)
            traps_g= traps_ar_g(i,j)
            tw_g= tw_ar_g(i,j)
            twcc_g= twcc_ar_g(i,j)
            mann_g= mann_ar_g(i,j)  !1.0/sk(i,j)
            manncc_g= manncc_ar_g(i,j) ! 1.0/skCC1(i,j)
            hbf_tz= (tw_g - bo_g)/(2.0*traps_g) !* bankfull depth

            if (tzeq_flag_g==0) then
            !* use x-sec attribute lookup tables
                icol= 2 !<-5
                co(i)= q_sk_multi*intp_elev(icol, i, j, xt)
                icol= 3 !<-6
                width= intp_elev(icol, i, j, xt)
            elseif (tzeq_flag_g==1) then
            !* use trapezoidal x-sec equations
                depth_tz= xt - z_g !z(i,j)
                call areacalc(depth_tz, area_tz)
                call hydRcalc(depth_tz, area_tz, hydr_tz)
                call Kcalc(depth_tz, area_tz, hydr_tz, conv_tz)
                if (depth_tz>hbf_tz) then
                    topwd_tz= twcc_g
                else
                    topwd_tz= bo_g + 2.0*traps_g*depth_tz
                end if
                co(i)= q_sk_multi*conv_tz
                width=topwd_tz
            endif

            sfi = ( q_g(i,j) / co(i) ) ** 2.0
            ybf_tz= hbf_tz + z_g !z_ar_g(i,j) !* bankfull depth/elevation for trapezoidal main channel
            if (xt.le.ybf_tz) then
            !* use Manning's N for main channel
                celty_g(i,j)= 5.0 / 3.0 * sfi ** 0.3 * abs(q_g(i,j)) ** 0.4 / width ** 0.4 / ( mann_g*(1/q_sk_multi)) ** 0.6
                emann_tz= -100.0 !* test only
            else
            !* when overbankflow, use equivalent manning's N
                !dsc= q_g(ts+1,i,j)
                depth_tz= xt - z_g !z_ar_g(i,j)
                emann_tz= emann_tmrf(depth_tz)
                celty_g(i,j)=5.0 / 3.0 * sfi ** 0.3 * abs(q_g(i,j)) ** 0.4 / width ** 0.4 / ( emann_tz*(1/q_sk_multi)) ** 0.6
            end if
            sum1= sum1 + celty_g(i,j)
            diffty_g(i,j) = abs(q_g(i,j)) / 2.0 / width / sfi
            sum2= sum2+ diffty_g(i,j)

            if (y_opt_g==2) then
            !* run diffusive depth computation
                elv_g(i,j) = elv_g(i+1,j) + sign ( sfi, q_g(i,j) ) * dx_ar_g(i,j)
            endif


!                write(301,"(I10,I2,I10,10F15.3)") ts, i, j, q_g(i,j), xt, depth_tz, co(i), sfi, width,&
!                                                emann_tz, celty_g(i,j), diffty_g(i,j), elv_g(i,j)

        enddo   !*do i=ncomp,2,-1

        !* averaged celty_g thorough all segments to prevent a wavefront of a cell with too high speed which
        !* produces wave breaking and instability
                !write(301,"(A10,A2,A10,2A15)") "ts", "i", "j", "celty_g(i,j)", "diffty_g(i,j)"
        celav= sum1/real(ncomp-1,KIND(sum1))
        diffav= sum2/real(ncomp-1,KIND(sum2))
        do i=1,ncomp
            celty_g(i,j) = celav
            diffty_g(i,j)= diffav
            if (diffty_g(i,j)>10.0) diffty_g(i,j)=10.0
                !write(301,"(I10,I2,I10,2F15.3)") ts, i, j, celty_g(i,j), diffty_g(i,j)
        end do
        cel_av_g(j)=celav


        deallocate(co)
        deallocate(harr_m, qarr_m, harr_f, qarr_f)
        deallocate(elvarr,  qarr)

    end subroutine elv_calc

!*----------------------------------------------------
!           Interpolate channel property by elev
!
!*----------------------------------------------------
    real(KIND=dp) function intp_elev(icol, inode, j, x)
        implicit none
        integer(KIND=i4b), intent(in) :: icol, inode, j
        real(KIND=dp), intent(in) :: x
        real(KIND=dp), allocatable :: dmyv(:)
        integer(KIND=i4b) ::  i1, irow, nrow
        real(KIND=dp) :: x1, y1, x2, y2, y

        nrow= nel_g
        allocate(dmyv(nrow))
        do i1=1, nrow
            dmyv(i1)= xsec_attr_rch_g(inode, j, i1, 1 )  !* elevation vector
        end do
        irow= locate(dmyv, x)
        if (irow.eq.0) irow= 1
        if (irow.eq.nrow) irow= nrow-1
        x1= xsec_attr_rch_g(inode, j, irow, 1 )
        y1= xsec_attr_rch_g(inode, j, irow, icol)
        x2= xsec_attr_rch_g(inode, j, irow+1, 1)
        y2= xsec_attr_rch_g(inode, j, irow+1, icol)

        y= LInterpol(x1,y1,x2,y2,x)
        intp_elev= y
        deallocate(dmyv)

    end function intp_elev

end module diff

