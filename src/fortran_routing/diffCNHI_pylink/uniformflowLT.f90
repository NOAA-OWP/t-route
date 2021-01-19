module flowLT

    !use arrays
    use var
    use subtools
    use nrtype

    implicit none

contains
!*-----------------------------------------------------------------------
!*          Uniform flow lookup table for
!*                           trap.main and rec.floodplain of a given node

!*  Create two types of lookup tables for a given node, p138-9,RM4:
!*  Table 1: normal depth vs. uniform flow in main channel
!*  Table 2: normal depth vs. uniform flow in flood plain

!*  Output -> ufhlt_m_g & ufqlt_m_g for Table 1; ufhlt_f_g & ufqlt_f_g for Table 2
!*------------------------------------------------------------------------
    subroutine uniflowLT_tz_alrch
    !subroutine uniflowLT_tz(jrch)
    !subroutine uniflowLT_tz(inode, jrch)
        implicit none

        !integer(KIND=i4b), intent(in) :: inode, jrch
        integer(KIND=i4b) :: i, i1, ncomp,j
        real(KIND=dp) :: hbf_tz, hincr_m, hincr_f, hstep, ufQ, hmx

        open(unit=401, file="./output/unifom flow LT m.txt", status="unknown")
        open(unit=402, file="./output/unifom flow LT f.txt", status="unknown")
        write(401,"(4A7, 2A15)") "ufLT_m","i","j","i1", "uf dsc_m", "nm depth_m"
        write(402,"(4A7, 2A15)") "ufLT_f","i","j","i1", "uf dsc_f", "nm depth_f"

        do j=1, nrch_g
            ncomp=frnw_g(j,1)
            do i=1, ncomp
                !if (tzeq_flag==1) then
                    !z0= z_ar_g(i) !z(i,j)
!                    bo0= bo_r(i) !bo(i,j)
!                    traps0= traps_r(i) !traps(i,j)
!                    tw0= tw_r(i) !Tw(i,j)
!                    twcc0= twcc_r(i) !TwCC(i,j)
!                    mann= mann_r(i) !mna(i,j)  !1.0/sk(i,j)
!                    manncc= manncc_r(i)   !mncca(i,j) ! 1.0/skCC1(i,j)
!                    So0= So_r(i)

                    !z0= z_ar_g(i,j)
                    bo_g= bo_ar_g(i,j)
                    traps_g= traps_ar_g(i,j)
                    tw_g= tw_ar_g(i,j)
                    twcc_g= twcc_ar_g(i,j)
                    mann_g= mann_ar_g(i,j)
                    manncc_g= manncc_ar_g(i,j)
                    So_g= So_ar_g(i,j)

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
                        call ufQ_tmrf(hstep, ufQ)
                        !* normal depth
                        ufhlt_m_g(i, j, i1)= hstep !ufhlt_m_g(i1, inode, jrch)= hstep
                        !* uniform discharge
                        ufqlt_m_g(i, j, i1)= ufQ !ufqlt_m_g(i1, inode, jrch)= ufQ

                        hstep= hstep + hincr_m
                        i1= i1+1
                    end do
                    hstep= hbf_tz
                    call ufQ_tmrf(hstep, ufQ)
                    i1= nhincr_m_g
                    !* normal depth
                    ufhlt_m_g(i, j, i1)= hstep !ufhlt_m_g(i1, inode, jrch)= hstep
                    !* uniform discharge
                    ufqlt_m_g(i, j, i1)= ufQ !ufqlt_m_g(i1, inode, jrch)= ufQ
                        !*test
                        do i1=1, nhincr_m_g
                            write(401,"(A7, 3I7, 2f15.3)") "ufLT_m", i, j, i1, ufqlt_m_g(i, j, i1), ufhlt_m_g(i, j, i1)
                        enddo
                        write(401,*)
                !*-------------------------------------------
                !* Lookup table for rectangular floodplain

                !*-------------------------------------------
                    hmx= timesDepth_g*hbf_tz
                    hincr_f= (hmx-hbf_tz)/real(nhincr_f_g, KIND(hbf_tz))
                    i1=1
                    hstep= hbf_tz
                    do while ((hstep < hmx).and.(i1<nhincr_f_g))
                    !do while (hstep < hmx)
                        call ufQ_tmrf(hstep, ufQ)
                        !* normal depth
                        ufhlt_f_g(i, j, i1)= hstep !ufhlt_f_g(i1, inode, jrch)= hstep
                        !* uniform discharge
                        ufqlt_f_g(i, j, i1)= ufQ !ufqlt_f_g(i1, inode, jrch)= ufQ

                        hstep= hstep + hincr_f
                        i1= i1+1
                    end do
                    hstep= hmx
                    call ufQ_tmrf(hstep, ufQ)
                    i1= nhincr_f_g
                    !* normal depth
                    ufhlt_f_g(i, j, i1)= hstep !ufhlt_f_g(i1, inode, jrch)= hstep
                    !* uniform discharge
                    ufqlt_f_g(i, j, i1)= ufQ !ufqlt_f_g(i1, inode, jrch)= ufQ
                         !*test
                        do i1=1, nhincr_f_g
                            write(402,"(A7, 3I7, 2f15.3)") "ufLT_f", i, j, i1, ufqlt_f_g(i, j, i1), ufhlt_f_g(i, j, i1)
                        enddo
                        write(402,*)
                !endif !* if (tzeq_flag==1) then
            enddo !*do i=1, ncomp
        enddo !*do j=1, nrch

    end subroutine uniflowLT_tz_alrch
!*-----------------------------------------------------------------------
!*          Uniform flow lookup table for
!*                           trap.main and rec.floodplain of a given node

!*  Create two types of lookup tables for a given node, p138-9,RM4:
!*  Table 1: normal depth vs. uniform flow in main channel
!*  Table 2: normal depth vs. uniform flow in flood plain

!*  Output -> ufhlt_m_g & ufqlt_m_g for Table 1; ufhlt_f_g & ufqlt_f_g for Table 2
!*------------------------------------------------------------------------
!    subroutine uniflowLT_tz
!    !subroutine uniflowLT_tz(jrch)
!    !subroutine uniflowLT_tz(inode, jrch)
!        implicit none
!
!        !integer(KIND=i4b), intent(in) :: inode, jrch
!        integer(KIND=i4b) :: i, i1
!        real(KIND=dp) :: hbf_tz, hincr_m, hincr_f, hstep, ufQ, hmx
!
!        open(unit=401, file="./output/unifom flow LT 1.txt", status="unknown")
!        open(unit=402, file="./output/unifom flow LT 2.txt", status="unknown")
!        write(401,"(2A7, 2A15)") "ufLT1","i","uf dsc_m", "nm depth_m"
!        write(402,"(2A7, 2A15)") "ufLT2","i","uf dsc_f", "nm depth_f"
!
!        do i=1, ncomp_unif
!            !if (tzeq_flag==1) then
!                !z0= z(inode, jrch)
!                !bo0= bo(inode, jrch)
!                !traps0= traps(inode, jrch)
!                !tw0= Tw(inode, jrch)
!                !twcc0= TwCC(inode, jrch)
!                !mann= mna(inode, jrch)
!                !manncc= mncca(inode, jrch
!                !So0= chbtslp(inode, jrch)
!                !* bankfull depth
!                !hbf_tz= (tw0 - bo0)/(2.0*traps0)
!                !z0= z_ar_g(i) !z(i,j)
!                bo0= bo_r(i) !bo(i,j)
!                traps0= traps_r(i) !traps(i,j)
!                tw0= tw_r(i) !Tw(i,j)
!                twcc0= twcc_r(i) !TwCC(i,j)
!                mann= mann_r(i) !mna(i,j)  !1.0/sk(i,j)
!                manncc= manncc_r(i)   !mncca(i,j) ! 1.0/skCC1(i,j)
!                So0= So_r(i)
!                hbf_tz= (tw0 - bo0)/(2.0*traps0) !* bankfull depth
!            !*-------------------------------------------
!            !* Lookup table for trapezoidal main channel
!
!            !*-------------------------------------------
!                !* depth interval
!                hincr_m= hbf_tz/real(nhincr_m_g, KIND(hbf_tz))
!                i1=1
!                hstep=0.0
!                do while (hstep < hbf_tz)
!                    call ufQ_tmrf(hstep, ufQ)
!                    !* normal depth
!                    ufhlt_mr(i1, i)= hstep !ufhlt_m_g(i1, inode, jrch)= hstep
!                    !* uniform discharge
!                    ufqlt_mr(i1, i)= ufQ !ufqlt_m_g(i1, inode, jrch)= ufQ
!
!                    hstep= hstep + hincr_m
!                    i1= i1+1
!                end do
!                hstep= hbf_tz
!                call ufQ_tmrf(hstep, ufQ)
!                i1= nhincr_m_g
!                !* normal depth
!                ufhlt_mr(i1, i)= hstep !ufhlt_m_g(i1, inode, jrch)= hstep
!                !* uniform discharge
!                ufqlt_mr(i1, i)= ufQ !ufqlt_m_g(i1, inode, jrch)= ufQ
!                    !*test
!                    do i1=1, nhincr_m_g
!                    write(401,"(A7, I7, 2f15.3)") "ufLT1", i, ufqlt_mr(i1, i), ufhlt_mr(i1, i)
!                    enddo
!                    !write(401,*)
!            !*-------------------------------------------
!            !* Lookup table for rectangular floodplain
!
!            !*-------------------------------------------
!                hmx= timesDepth*hbf_tz
!                hincr_f= (hmx-hbf_tz)/real(nhincr_f_g, KIND(hbf_tz))
!                i1=1
!                hstep= hbf_tz
!                do while (hstep < hmx)
!                    call ufQ_tmrf(hstep, ufQ)
!                    !* normal depth
!                    ufhlt_fr(i1, i)= hstep !ufhlt_f_g(i1, inode, jrch)= hstep
!                    !* uniform discharge
!                    ufqlt_fr(i1, i)= ufQ !ufqlt_f_g(i1, inode, jrch)= ufQ
!
!                    hstep= hstep + hincr_f
!                    i1= i1+1
!                end do
!                hstep= hmx
!                call ufQ_tmrf(hstep, ufQ)
!                i1= nhincr_f_g
!                !* normal depth
!                ufhlt_fr(i1, i)= hstep !ufhlt_f_g(i1, inode, jrch)= hstep
!                !* uniform discharge
!                ufqlt_fr(i1, i)= ufQ !ufqlt_f_g(i1, inode, jrch)= ufQ
!                     !*test
!                    do i1=1, nhincr_f_g
!                    write(402,"(A7, I7, 2f15.3)") "ufLT2", i, ufqlt_fr(i1, i), ufhlt_fr(i1, i)
!                    enddo
!!                    write(402,*)
!            !endif !* if (tzeq_flag==1) then
!        enddo !*do i=1, ncomp_unif
!
!    end subroutine uniflowLT_tz

end module flowLT
