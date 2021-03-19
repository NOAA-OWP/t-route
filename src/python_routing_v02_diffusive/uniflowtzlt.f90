module flowlt
    implicit none
    double precision :: bo_g, traps_g, tw_g, twcc_g, mann_g, manncc_g, So_g

contains
!*-----------------------------------------------------------------------
!*          Uniform flow lookup table for
!*                           trap.main and rec.floodplain of a given node

!*  Create two types of lookup tables for a given node, p138-9,RM4:
!*  Table 1: normal depth vs. uniform flow in main channel
!*  Table 2: normal depth vs. uniform flow in flood plain

!*  Output -> ufhlt_m_g & ufqlt_m_g for Table 1; ufhlt_f_g & ufqlt_f_g for Table 2
!*------------------------------------------------------------------------
    subroutine uniflowtzlt(mxncomp_g, nrch_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g, &
                                    mann_ar_g, manncc_ar_g, so_ar_g, &
                                    nhincr_m_g, nhincr_f_g, ufhlt_m_g, ufqlt_m_g, ufhlt_f_g, ufqlt_f_g,&
                                    frnw_col, dfrnw_g, timesdepth_g)
        implicit none
        integer, intent(in) :: mxncomp_g, nrch_g
        integer, intent(in) :: nhincr_m_g, nhincr_f_g, frnw_col
        double precision, dimension(mxncomp_g, nrch_g), intent(in) :: bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g
        double precision, dimension(mxncomp_g, nrch_g), intent(in) :: mann_ar_g, manncc_ar_g, so_ar_g
        double precision, dimension(nrch_g, frnw_col), intent(in) :: dfrnw_g 
        double precision, intent(in) :: timesdepth_g
        double precision, dimension(mxncomp_g, nrch_g, nhincr_m_g), intent(out) :: ufhlt_m_g,  ufqlt_m_g
        double precision, dimension(mxncomp_g, nrch_g, nhincr_f_g), intent(out) :: ufhlt_f_g, ufqlt_f_g
        integer :: i, i1, ncomp,j
        double precision :: hbf_tz, hincr_m, hincr_f, hstep, ufQ, hmx, incr_m
        integer, dimension(:,:), allocatable :: frnw_g
        allocate(frnw_g(nrch_g, frnw_col))
        frnw_g=int(dfrnw_g) 

        !open(unit=401, file="./output/unifom flow LT m.txt", status="unknown")
        !open(unit=402, file="./output/unifom flow LT f.txt", status="unknown")
        !write(401,"(4A7, 2A15)") "ufLT_m","i","j","i1", "uf dsc_m", "nm depth_m"
        !write(402,"(4A7, 2A15)") "ufLT_f","i","j","i1", "uf dsc_f", "nm depth_f"

        do j=1, nrch_g
            ncomp=frnw_g(j,1)
            do i=1, ncomp
                bo_g= bo_ar_g(i,j)
                traps_g= traps_ar_g(i,j)
                tw_g= tw_ar_g(i,j)
                twcc_g= twcc_ar_g(i,j)
                mann_g= mann_ar_g(i,j)
                manncc_g= manncc_ar_g(i,j)
                So_g= so_ar_g(i,j)

                hbf_tz= (tw_g - bo_g)/(2.0*traps_g) !* bankfull depth
                !*-------------------------------------------
                !* Lookup table for trapezoidal main channel

                !*-------------------------------------------
                !* depth interval
                hincr_m= hbf_tz/real(nhincr_m_g, KIND(hbf_tz))
                i1=1
                hstep=0.0
                do while ((hstep < hbf_tz).and.(i1<nhincr_m_g))
                !do while (hstep < hbf_tz)
                    call ufq_tmrf(hstep, ufQ)
                    !* normal depth
                    ufhlt_m_g(i, j, i1)= hstep 
                    !* uniform discharge
                    ufqlt_m_g(i, j, i1)= ufQ 

                    hstep= hstep + hincr_m
                    i1= i1+1
                end do
                hstep= hbf_tz
                call ufq_tmrf(hstep, ufQ)
                i1= nhincr_m_g
                !* normal depth
                ufhlt_m_g(i, j, i1)= hstep 
                !* uniform discharge
                ufqlt_m_g(i, j, i1)= ufQ 
                        !*test
                        !do i1=1, nhincr_m_g
                        !    write(401,"(A7, 3I7, 2f15.3)") "ufLT_m", i, j, i1, ufqlt_m_g(i, j, i1), ufhlt_m_g(i, j, i1)
                        !enddo
                        !write(401,*)
                !*-------------------------------------------
                !* Lookup table for rectangular floodplain

                !*-------------------------------------------
                hmx= timesdepth_g*hbf_tz
                hincr_f= (hmx-hbf_tz)/real(nhincr_f_g, KIND(hbf_tz))
                i1=1
                hstep= hbf_tz
                do while ((hstep < hmx).and.(i1<nhincr_f_g))
                    call ufq_tmrf(hstep, ufQ)
                    !* normal depth
                    ufhlt_f_g(i, j, i1)= hstep 
                    !* uniform discharge
                    ufqlt_f_g(i, j, i1)= ufQ 
                    hstep= hstep + hincr_f
                    i1= i1+1
                 end do
                 hstep= hmx
                 call ufq_tmrf(hstep, ufQ)
                 i1= nhincr_f_g
                 !* normal depth
                 ufhlt_f_g(i, j, i1)= hstep 
                 !* uniform discharge
                 ufqlt_f_g(i, j, i1)= ufQ 
                         !*test
                        !do i1=1, nhincr_f_g
                        !    write(402,"(A7, 3I7, 2f15.3)") "ufLT_f", i, j, i1, ufqlt_f_g(i, j, i1), ufhlt_f_g(i, j, i1)
                        !enddo
                        !write(402,*)
            enddo !*do i=1, ncomp
        enddo !*do j=1, nrch
        deallocate(frnw_g) 
    end subroutine uniflowtzlt
!*----------------------------------------------------------------------
!*          Uniform flow discharge for trapz. main and rect. floodplain

!*  p.102, RM4
!*----------------------------------------------------------------------
    subroutine ufq_tmrf(hxs, ufQ)
        implicit none

        double precision, intent(in) :: hxs
        double precision, intent(out) :: ufQ
        double precision :: Axs, hydR, WP, Twd, TCCwd, So
        double precision :: bwd, sslp, hbf, ufQ1, ufQ2, ufQ3

        bwd= bo_g
        sslp= traps_g
        Twd= tw_g       !* top width of main channel
        TCCwd= twcc_g   !* top width of compound channel
        So= So_g

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

    end subroutine ufq_tmrf
end module flowlt
