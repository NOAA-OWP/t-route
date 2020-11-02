!Written by Dong Ha Kim, NOAA's Office of Water Prediction, National Water Center
module flowLT

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

!*  Output -> ufhlt_m & ufqlt_m for Table 1; ufhlt_f & ufqlt_f for Table 2
!*------------------------------------------------------------------------
    subroutine uniflowLT_tz

        implicit none

        integer(KIND=i4b) :: i, i1
        real(KIND=dp) :: hbf_tz, hincr_m, hincr_f, hstep, ufQ, hmx

        !open(unit=401, file="./output/unifom flow LT 1.txt", status="unknown")
        !open(unit=402, file="./output/unifom flow LT 2.txt", status="unknown")

        do i=1, ncomp
            bo0= bo_ar(i)
            traps0= traps_ar(i)
            tw0= tw_ar(i)
            twcc0= twcc_ar(i)
            mann= mann_ar(i)
            manncc= manncc_ar(i)
            So0= So_ar(i)
            hbf_tz= (tw0 - bo0)/(2.0*traps0) !* bankfull depth
        !*-------------------------------------------
        !* Lookup table for trapezoidal main channel

        !*-------------------------------------------
            !* depth interval
            hincr_m= hbf_tz/real(nhincr_m, KIND(hbf_tz))
            i1=1
            hstep=0.0
            do while (hstep < hbf_tz)
                call ufQ_tmrf(hstep, ufQ)
                !* normal depth
                ufhlt_m(i1, i)= hstep
                !* uniform discharge
                ufqlt_m(i1, i)= ufQ

                hstep= hstep + hincr_m
                i1= i1+1
            end do
            hstep= hbf_tz
            call ufQ_tmrf(hstep, ufQ)
            i1= nhincr_m
            !* normal depth
            ufhlt_m(i1, i)= hstep
            !* uniform discharge
            ufqlt_m(i1, i)= ufQ
        !*-------------------------------------------
        !* Lookup table for rectangular floodplain

        !*-------------------------------------------
            hmx= timesDepth*hbf_tz
            hincr_f= (hmx-hbf_tz)/real(nhincr_f, KIND(hbf_tz))
            i1=1
            hstep= hbf_tz
            do while (hstep < hmx)
                call ufQ_tmrf(hstep, ufQ)
                !* normal depth
                ufhlt_f(i1, i)= hstep
                !* uniform discharge
                ufqlt_f(i1, i)= ufQ

                hstep= hstep + hincr_f
                i1= i1+1
            end do
            hstep= hmx
            call ufQ_tmrf(hstep, ufQ)
            i1= nhincr_f
            !* normal depth
            ufhlt_f(i1, i)= hstep
            !* uniform discharge
            ufqlt_f(i1, i)= ufQ
        enddo !*do i=1, ncomp

    end subroutine uniflowLT_tz

end module flowLT
