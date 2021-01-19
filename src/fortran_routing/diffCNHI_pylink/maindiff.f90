!            FINITE DIFFERENCE METHOD
!
!  A program for one dimensional flow in open channel
!
!** MSH implementation based on p103~109,RM4
!* Input channel data == Pocono sampleRouteLink3

!** Changes against diffusive_irregular_CNHI_mtime_mrch **
!* From:
!    allocate(celty_g(ntss_g, mxncomp_g, nlinks), diffty_g(ntss_g, mxncomp_g, nlinks))
!    allocate(q_g(ntss_g, mxncomp_g, nlinks), qpx_g(ntss_g, mxncomp_g, nlinks))
!    allocate(elv_g(ntss_g, mxncomp_g, nlinks))

!* To:
!    allocate(celty_g(mxncomp_g, nlinks), diffty_g(mxncomp_g, nlinks))
!    allocate(q_g(mxncomp_g, nlinks), qpx_g(mxncomp_g, nlinks))
!    allocate(elv_g(mxncomp_g, nlinks))

!** Changes against diffusive_irregular_CNHI_mtime_mrch_slimvar_g **
!* -> Instead of varying dt, dx is modified at each simulation time step in order to make
!*    Courant number less than 1.  Refer to p.2, RM5.


!program main2
module mdiff

    use var
    use subtools
    use diff
    !use simtime
    use nrtype
    use flowLT
    use mdxtools

    implicit none



contains
    !*--------------------------------------------------------------------------------
    !*       Route network by using Python-to-Fortran network traversal map together
    !*          with diffusive routing engines
    !
    !*--------------------------------------------------------------------------------
    subroutine routenw

        ! Local storage
        integer(KIND=i4b) :: ts, n, i, j, ts_ev, i1
        integer(KIND=i4b) :: ncomp

        real(KIND=dp) ::  tc, tf0
        real(KIND=dp) :: qlatf


        integer(KIND=i4b) :: nusrch, rch, usrchj, ncomp_usrchj, dsrchj
        integer(KIND=i4b) :: lntss, lmxncomp, lnrch
        real(KIND=dp) :: qjt
        real(KIND=dp), dimension(:), allocatable ::tarr_ql, varr_ql, tarr_ub, varr_ub, tarr_db, varr_db

        integer(KIND=i4b) :: mi
        real(KIND=dp) :: vr, vr2, mtp, mtp2, mn, mn2, yintc, cel_mx, rmnd, saveInterval_min
        integer(KIND=i4b) :: dmyi, dmyi1, dmyi2 !test only
        real(KIND=dp) :: dmyr, dmyr1, dmyr2, dmyr3, dmyr4 !test only

        integer(KIND=i4b) :: nldxj, mj, ncomp_m, i_m, x_mflag, i_s, i2, xflag, i_m_s
        real(KIND=dp) ::  delxj_mn, delxj_m, sumdxj, x_m, sumx1, sumx2, dx_s, dx_i_s, delx, u1, u2, u_m
        real(KIND=dp) :: sumqlatdx, sumqlatx_m, sumqlatdx_m, qlatf_m, x_o, sumx1_m, sumx2_m, delx_m, u_o
        real(KIND=dp) :: so_us, c_us, so_ds, c_ds, adjso, initval
        integer(KIND=i4b) :: frj, iseg

!*Py        nrch_g=737
!*Py         mxncomp_g=7

!*Py         dtini_g= 300.0         !* initial time interval in [sec]
!*Py         t0_g= 0.0             !* simulation starting time in [hr]
!*Py         tfin_g= 2.0       !* simulation ending time in [hr] #3670.0
!*Py         saveInterval_g=dtini_g    !* output publishing time interval in [sec]
!*Py         cfl_g= 0.99            !* max. allowable Courant number
!*Py         dt_ql_g= 3600.0  !* instead of dtinput #time interval of qlateral input data from wrf-hydro [sec]
!*Py         dt_ub_g= dt_ql_g  !* time interval of input data for upstream boundary condition  [sec]
!*Py         dt_db_g=  900.0 !* time interval of input data for downstream boundary condition  [sec]
!*Py         saveInterval_ev_g=saveInterval_g !for evaluation later

!*Py         nl_ubcd_ar_g=1116 !* row number of ubcd_ar.txt

!*Py         ntss_ev_g= int((tfin_g - t0_g)*3600.0/saveInterval_ev_g, KIND(ntss_ev_g))+4
!*Py         ntss_g= 40*ntss_ev_g  !* int(10.0*(tfin - t0)*3600.0/dtini, KIND(ntss_g))

!*Py         nts_ql_g= int( (tfin_g - t0_g)*3600.0/dt_ql_g+1, KIND(nts_ql_g) ) !* based on time of [sec]
!*Py         nts_ub_g= nts_ql_g
!*Py         nts_db_g= int( (tfin_g - t0_g)*3600.0/dt_db_g+1, KIND(nts_db_g) ) !* based on time of [sec]

        allocate(celty_g(mxncomp_g, nrch_g), diffty_g(mxncomp_g, nrch_g))
        allocate(q_g(mxncomp_g, nrch_g), elv_g(mxncomp_g, nrch_g), qpx_g(mxncomp_g, nrch_g))

!*Py         allocate(qlat_g(nts_ql_g,mxncomp_g,nrch_g))
!*Py         allocate(q_ev_g(ntss_ev_g, mxncomp_g, nrch_g), elv_ev_g(ntss_ev_g, mxncomp_g, nrch_g))

        allocate(tarr_ql(nts_ql_g), varr_ql(nts_ql_g))
        allocate(tarr_ub(nts_ub_g), varr_ub(nts_ub_g))
        allocate(tarr_db(nts_db_g), varr_db(nts_db_g))
        allocate(cel_av_g(nrch_g))

!*Py         nhincr_m_g=10
!*Py         nhincr_f_g=20
!*Py         allocate(ufhlt_m_g(mxncomp_g,nrch_g,nhincr_m_g), ufqlt_m_g(mxncomp_g,nrch_g,nhincr_m_g))
!*Py         allocate(ufhlt_f_g(mxncomp_g,nrch_g,nhincr_f_g), ufqlt_f_g(mxncomp_g,nrch_g,nhincr_f_g))



        ! Allocate arrays
        !call setup_arrays(ntss_g, ntsi, mxncomp_g, nlinks, mxnbrch, nends)
        open(unit=1, file="./output/q_ev_g elv_ev_g.txt")
        open(unit=401, file="./output/f2py_FallsLake_diffusive_irregular_CNHI_mtime&
                        _mrch_slimvar_g_mdfdx.txt")
!        open(unit=2, file="./output/q elv.txt")
!        open(unit=11, file="./output/frnw_g.txt")
        open(unit=12, file="./output/adjust so.txt")
!        open(unit=13, file="./output/adj dx test.txt")
!        open(unit=91, file="./output/diagonastic.txt", status="unknown")
!        open(unit=311, file="./output/validate qlat_ar import.txt")
!        open(unit=321, file="./output/validate ubcd_ar import.txt")
!
!
!        open(unit=21, file="./input/frnw_ar.txt")
!        open(unit=22, file="./input/z_bo_traps_tw_twcc_m_mcc_so_dx.txt")
!        open(unit=22, file="./input/z_ar.txt")
!        open(unit=23, file="./input/bo_ar.txt")
!        open(unit=24, file="./input/traps_ar.txt")
!        open(unit=25, file="./input/tw_ar.txt")
!        open(unit=26, file="./input/twcc_ar.txt")
!        open(unit=27, file="./input/mann_ar.txt")
!        open(unit=28, file="./input/manncc_ar.txt")
!        open(unit=29, file="./input/so_ar.txt")
!        open(unit=30, file="./input/dx_ar.txt")
!        open(unit=31, file="./input/qlat_ar.txt")
!        open(unit=32, file="./input/ubcd_ar.txt")
!        open(unit=33, file="./input/dbcd_ar.txt")

!        allocate(frnw_g(nrch_g,8))
!        do j=1, nrch_g
!            read(21,*) (frnw_g(j,i),i=1,8)
!        enddo
!
!        write(11,"(6A15)") "j","frnw_g(j,1)", "frnw_g(j,2)", "frnw_g(j,3)", "frnw_g(j,4)", "frnw_g(j,5)"
!        do j=1, nrch_g
!            write(11,"(6I15)") j,frnw_g(j,1),frnw_g(j,2),frnw_g(j,3),frnw_g(j,4),frnw_g(j,5)
!        enddo


!        allocate( z_ar_g(mxncomp_g, nrch_g), bo_ar_g(mxncomp_g, nrch_g), traps_ar_g(mxncomp_g, nrch_g) )
!        allocate( tw_ar_g(mxncomp_g, nrch_g), twcc_ar_g(mxncomp_g, nrch_g) )
!        allocate( mann_ar_g(mxncomp_g, nrch_g), manncc_ar_g(mxncomp_g, nrch_g) )
!        allocate( So_ar_g(mxncomp_g, nrch_g), dx_ar_g(mxncomp_g, nrch_g) )
!
!
!        do i=1, mxncomp_g
!            read(22,*) (z_ar_g(i,j),j=1,nrch_g)
!        enddo
!        do i=1, mxncomp_g
!            read(23,*) (bo_ar_g(i,j),j=1,nrch_g)
!        enddo
!        do i=1, mxncomp_g
!            read(24,*) (traps_ar_g(i,j),j=1,nrch_g)
!        enddo
!        do i=1, mxncomp_g
!            read(25,*) (tw_ar_g(i,j),j=1,nrch_g)
!        enddo
!        do i=1, mxncomp_g
!            read(26,*) (twcc_ar_g(i,j),j=1,nrch_g)
!        enddo
!        do i=1, mxncomp_g
!            read(27,*) (mann_ar_g(i,j),j=1,nrch_g)
!        enddo
!        do i=1, mxncomp_g
!            read(28,*) (manncc_ar_g(i,j),j=1,nrch_g)
!        enddo
!        do i=1, mxncomp_g
!            read(29,*) (So_ar_g(i,j),j=1,nrch_g)
!        enddo
!        do i=1, mxncomp_g
!            read(30,*) (dx_ar_g(i,j),j=1,nrch_g)
!        enddo



        !++-----------------------------------------------------------------------------------+
        !+        adjust abnormally small channel bottom slope
        !+ Definition of abnormal slope:
        !+  1. slopes less than a chosen lower limit (=so_llm)
        !+ Repair method:
        !+  1. take average of slopes of adjacent segments
        !++-----------------------------------------------------------------------------------+
        write(12,"(4A15)") "bf_adj", "i", "j", "So_ar_g(i,j)"
        do j=1,nrch_g
            ncomp=frnw_g(j,1)
            do i=1, ncomp-1
                write(12,"(A15, 2I15, F15.6)") "bf_adj", i, j, So_ar_g(i,j)
            enddo
        enddo
!        write(12,*)

!*Py         so_llm_g=0.0005
        do j=1,nrch_g
            ncomp=frnw_g(j,1)
            do i=1, ncomp-1
                if (So_ar_g(i,j).lt.so_llm_g) then
                    !* adjacent upstream segment's slope
                    so_us=0.0
                    c_us=0.0
                    i2=i-1
                    do while (i2.ge.1)
                        if (So_ar_g(i2,j).ge.so_llm_g) then
                            so_us= So_ar_g(i2,j)
                            c_us=1.0
                            exit
                        endif
                        i2=i2-1
                    end do
                    !* adjacent downstream segment's slope
                    so_ds=0.0
                    c_ds=0.0
                    i2=i+1
                    do while (i2.le.ncomp-1)
                        if (So_ar_g(i2,j).ge.so_llm_g) then
                            so_ds= So_ar_g(i2,j)
                            c_ds=1.0
                            exit
                        endif
                        i2=i2+1
                    end do
                    if (c_us+c_ds.gt.0.0) then
                        adjso= (so_us+so_ds)/(c_us+c_ds)
                        So_ar_g(i,j)= adjso
                    else
                        So_ar_g(i,j)= so_llm_g
                    endif
                endif
            enddo
        enddo

        write(12,"(4A15)") "af_adj", "i", "j", "So_ar_g(i,j)"
        do j=1,nrch_g
            ncomp=frnw_g(j,1)
            do i=1, ncomp-1
                write(12,"(A15, 2I15, F15.6)") "af_adj", i, j, So_ar_g(i,j)
            enddo
        enddo



        !++-----------------------------------------------------------------------------------+
        !+        uniform flow lookup tables

        !+  Output:
        !+    Lookup table 1 for main ch.: [ufhlt_m_g(i1,i,j), ufqlt_m_g(i1,i,j)] as
        !+                                  normal depth & uniform discharge for i1=1, nhincr_m_g
        !*    Lookup table 2 for floodplain: [ufhlt_f_g(i1,i,j), ufqlt_f_g(i1,i,j)] as
        !+                                  normal depth & uniform discharge for i1=1, nhincr_f_g
        !++-----------------------------------------------------------------------------------+
!*Py         timesDepth_g=5.0
!*Py         call uniflowLT_tz_alrch


    !    open(unit=41, file="./output/qlat_g.txt", status="unknown")
    !    open(unit=42, file="./output/q_headsegment.txt", status="unknown")
    !    open(unit=43, file="./output/y_downstream_bd.txt", status="unknown")
    !    write(41,"(3A7,2A12)") "n","i","j","[mim]",  "qlat_g(n,i,j)"
        !++-------------------------------------------+
        !+              Lateral inflow
        !+
        !++-------------------------------------------
!        vr=2.5
!        mtp=5.0
!        do n=1,nts_ql_g
!            do j=1,nrch_g
!                ncomp=frnw_g(j,1)
!                mn= 40.0+ real(j,KIND(mn))
!                do i=1, ncomp-1
!                    !* qlat_g in m2/sec
!                    qlat_g(n,i,j)= mtp*exp(-0.5*((real(n,KIND(mn))-mn)/vr)**2.0)/(vr*(2.0*3.14)**0.5) + 0.1
!                    qlat_g(n,i,j)= qlat_g(n,i,j)/dx_ar_g(i,j)
!                    !write(41,"(3I7,2F12.7)") n,i,j,tarri(n), qlat_g(n,i,j)
!                end do
!                qlat_g(n,ncomp,j)= 0.0
!            enddo
!        enddo
!                write(311,"(3A15,A15, A25)") "j", "i", "n","tarr_ql(n)", "qlat_g(n,i,j)"
!        do j=1, nrch_g
!            ncomp=frnw_g(j,1)
!            do i=1, ncomp
!                do n=1, nts_ql_g
!                    read(31,*) frj, iseg, tarr_ql(n), qlat_g(n,i,j) !* qlat_g(n,i,j) in [m^2/sec]
!                    write(311,"(3I15,f15.2,f25.20)") j, i, n, tarr_ql(n), qlat_g(n,i,j)
!                end do
!            end do
!        enddo
        do n=1, nts_ql_g
            tarr_ql(n)= t0_g*60.0 + dt_ql_g*real(n-1,KIND(dt_ql_g))/60.0 !* [min]
        end do

        !++--------------------------------------------------------------+
        !+              Upstream boundary condition
        !+
        !+ Assume that hydrographs are always available at the upper ends
        !+ of the most upstream links.
        !++--------------------------------------------------------------+
        !write(42,"(A7,3A15)") "n","[min]", "q_g(n,1,1)","q_g(n,1,2)" !,"q_g(n,1,4)"
!*Py        allocate( ubcd_g(nts_ub_g, nrch_g) )
!        vr=3.0
!        mtp=20.0
!        mn= 10.0
!        yintc= 10.0
!        do j=1, nrch_g
!            if (frnw_g(j,3)==0) then !* frnw_g(j,3) indicates the number of upstream reaches.
!                yintc= yintc+0.2
!                do n=1,nts_ub_g
!                    ubcd_g(n,j)= mtp*exp(-0.5*((real(n,KIND(mn))-mn)/vr)**2.0)/(vr*(2.0*3.14)**0.5)+yintc
!                enddo
!            endif
!        enddo

!            write(321,"(3A15,2A20)") "row", "j", "i", "tarr_ub(i)", "ubcd_g(i,frj)"
!        do j=1, nrch_g
!            do n=1,nts_ub_g
!                ubcd_g(n,j)=-1.0 !* initialize
!            end do
!        end do
!        !* ubcd_g in [m^3/sec]
!        do i1=1, nl_ubcd_ar_g
!            read(32,*) frj, i, tarr_ub(i), ubcd_g(i,frj)
!                         write(321,"(3I15,f20.2,f20.10)") i1, frj, i, tarr_ub(i), ubcd_g(i,frj)
!        enddo
        do n=1, nts_ub_g
            tarr_ub(n)= t0_g*60.0 + dt_ub_g*real(n-1,KIND(dt_ub_g))/60.0 !* [min]
        enddo
       !++----------------------------------------------------------------+
        !+          Stage data at the downstream end
        !+
        !+
        !+ when flow conditions of lower end of link j are known.
        !+ Assume that stage time-series data is always available at the
        !+ lower end of the most downstream link.
        !++----------------------------------------------------------------+
!*Py        allocate( dbcd_g(nts_db_g) )
!        mtp2=5.0
!        mn2=mn+10.0
!        vr2=2.0*vr
!        ncomp=frnw_g(nrch_g,1)
!        do n=1, nts_db_g
!            dbcd_g(n) = mtp2*exp(-0.5*((real(n,KIND(mn2))-mn2)/vr2)**2.0)/(vr2*(2.0*3.14)**0.5)+1.0 + z_ar_g(ncomp,nrch_g)
!                !write(*,"(3I7,3F15.4)") n,ncomp,j,tarri(n), y(n,ncomp,j), y(n,ncomp,j)-z(ncomp,j)
!        end do
!        do n=1, nts_db_g
!            read(33,*) tarr_db(n), dbcd_g(n) !*  dbcd_g(n)  in [meter] based on NAD88 meter datum that RouteLink uses
!        end do
        do n=1, nts_db_g
            tarr_db(n)= t0_g*60.0 + dt_db_g*real(n-1,KIND(dt_db_g))/60.0 !* [min]
        enddo

        !+++---------------------------------------------------------------------------------+
        !+                  INITIAL CONDITION of q
        !+
        !+++---------------------------------------------------------------------------------+
        !+ 1) Headbasin link: initialize q along all the nodes with the measured data at n=1
        !+ until the link meets a junction.
        ts=1
        do j=1, nrch_g
            ncomp=frnw_g(j,1)
        !* For the head node of a reach,
            if (frnw_g(j,3)==0) then !* frnw_g(j,3) indicates the number of upstream reaches.
            !* head water reach
               i=1
               q_g(i,j)= ubcd_g(ts,j)
            else
            !* at a junction
                qjt= 0.0
                nusrch= frnw_g(j,3) !* then number of upstream reaches
                do rch=1, nusrch
                    usrchj= frnw_g(j,3+rch) !* js corresponding to upstream reaches
                    ncomp_usrchj= frnw_g(usrchj,1)
                    qjt= qjt + q_g(ncomp_usrchj, usrchj)
                enddo
                i=1
                q_g(i,j)= qjt
            endif
        !* For the following nodes within a reach after the head node,
            do i=2, ncomp
                !* qlat_g(n,i,j) in unit of m2/sec
                q_g(i,j) = q_g(1,j) + qlat_g(ts,i-1,j)*dx_ar_g(i-1,j)
            enddo
        enddo

    !            !* Output initial condition at n=1
    !            open(unit=81, file="./output/elev.txt", status="unknown")
    !            open(unit=82, file="./output/q.txt", status="unknown")
    !            open(unit=83, file="./output/celerity and inflow.txt", status="unknown")

    !
    !            write(81,"(4A10, 3A15)") "ts     ","   tc[min]", "i", "j", "elv(ts,i,j)",  "z(i,j)", "depth(i,j)"
    !            write(82,"(4A10, 4A15)") "ts     ","   tc[min]", "i", "j", "q_g(ts,i,j)",  "qlatj_g(i)", "qlatj_g*dx_ar_g", "dx_ar_g"
    !            write(83,"(4A10, 3A15)") "ts     ","   tc[min]", "i", "j", "celerity",  "qlat_g[csm]"

        tc = t0_g*60.0     !!! t0 is in hour. tc is in minutes
        ts=1    !*simulation time step in parallel with tc.
        ts_ev=1 !* time step for outputting q and elv at evaluation time
            write(1,"(4A10, 3A20)") "tc[min]", "ts_ev",  "i", "j", "q_ev_g", "elv_ev_g", "depth"
            do j=1, nrch_g
                ncomp= frnw_g(j,1)
                do i=1, ncomp
                    q_ev_g(ts_ev, i, j)=q_g(i,j)
                    write(1,"(f10.3, 3I10, F20.4)") tc, ts_ev, i, j, q_ev_g(ts_ev,i,j)
                enddo
            enddo
        ts_ev=ts_ev+1
!*Py        theta_g=1.0
!*Py         y_opt_g=1 !* 1 for normal depth(kinematic); 2 for dept of diffusive wave.
!*Py         tzeq_flag_g=1
        cel_av_g=1.0 !* initial value
        initval= -1.0 !* initial value after allocating variables



        do while ( tc .lt. tfin_g*60.0)

            !call cal_dt(tc, cel_mx) !output-> dtini [sec]
            !*-------------------------------------------------------------------------
            !* Compute q and qpx_g at ts+1 from upstream reach to downstream reach while
            !* computing from downstream node to upstream node for each reach.
            !*-------------------------------------------------------------------------
            do j=1, nrch_g
                !* Compute delxj_mn, p.2-3-3,RM5
                delxj_mn= cel_av_g(j)*dtini_g/cfl_g

                ncomp= frnw_g(j,1) !ncomp=nx1(j)

                !* interpolate qlat_g at the current time tc
                allocate(qlatj_g(ncomp))
                qlatj_g=initval !* initialize elements' dummy values
                do i=1,ncomp-1
                    do n=1,nts_ql_g
                        !tarr_ql(n)= dt_ql_g*real(n,KIND(dt_ql_g))/60.0 !* [min]
                        varr_ql(n)= qlat_g(n,i,j) !* qlat_g(n,i,j) in unit of m2/sec
                    enddo
                    qlatj_g(i)= intp_y(nts_ql_g, tarr_ql, varr_ql, tc)
                enddo

                if (ts==1) then
                !* initial value at ts=1
                    do i=1, ncomp
                        qpx_g(i,j)=0.0
                        celty_g(i,j)=1.0
                        diffty_g(i,j)=10.0
                    enddo
                end if

                !* Check how many segments of reach j are larger than delxj_mn
                nldxj=0
                do i=1, ncomp-1
                    if (dx_ar_g(i,j).ge.delxj_mn) then
                        nldxj=nldxj+1
                    end if
                end do
                !*-------------------------------------------------------------------------------------------
                !* Compute eei_g, ffi_g, exi_g, fxi_g from upstream to downstream nodes for a single j reach
                !*-------------------------------------------------------------------------------------------
                if (nldxj==ncomp-1) then
                !* all the segments of reach j are larger than delxj_mn, so execute sub. ef_calc as it is.
                    allocate(eei_g(ncomp), ffi_g(ncomp), exi_g(ncomp), fxi_g(ncomp))
                    !* upstream terminal node boundary condition
                    eei_g(1) = 1.0
                    ffi_g(1) = 0.
                    exi_g(1) = 0.
                    fxi_g(1) = 0.

                    !* compute eei_g, ffi_g, exi_g, and fxi_g for updating q and qpx_g at ts+1 after this call.
                     call ef_calc(ncomp, j)

!                        write(91,"(4A10, 10A20)") "ts","CURR:tc_M", "i", "j", "eei_g", "ffi_g","exi_g","fxi_g",&
!                                        "q_g","qpx_g","qlatj_g(i)", "celty_g", "diffty_g", "dx_ar_g(i,j)"
!                        do i=1,ncomp
!                            write(91,"(I10, F10.1, 2I10, 10F20.8)") ts, tc, i, j, eei_g(i), ffi_g(i), exi_g(i), fxi_g(i),&
!                                                q_g(i,j), qpx_g(i,j), qlatj_g(i), celty_g(i,j), diffty_g(i,j),&
!                                                dx_ar_g(i,j)
!                        enddo
                    !*----------------------------------------------------------------------------------
                    !* compute q and qpx_g at ts+1 from downstream to upstream nodes for a single j reach
                    !*----------------------------------------------------------------------------------
                    do n=1,nts_ql_g
                        varr_ql(n)= qlat_g(n,ncomp-1,j) !* qlat_g(n,i,j) in unit of m2/sec
                    enddo
                    tf0= tc+dtini_g/60.0
                    qlatf= intp_y(nts_ql_g, tarr_ql, varr_ql, tf0)
                    q_g(ncomp,j)=q_g(ncomp-1,j) + qlatf*dx_ar_g(ncomp-1,j) !* qlatf at ts+1 in unit of m2/sec

                    qpx_g(ncomp,j)=0.0
                    do i=ncomp-1, 1, -1
                        q_g(i,j) = eei_g(i) * q_g(i+1,j) + ffi_g(i)    !* q_g(i) at time t+dtini/60 [min] (=n+1)
                        qpx_g(i,j)= exi_g(i) *qpx_g(i+1,j) + fxi_g(i)  !* Qx(i) at time t+dtini/60 [min] (=n+1)
                    enddo
                    deallocate(eei_g, ffi_g, exi_g, fxi_g)

!                        write(91,"(4A10, 3A20)") "ts+1","FUT:tc_M", "i", "j", "q_g","qpx_g","qlatf"
!                        do i=1,ncomp
!                            write(91,"(I10, F10.1, 2I10, 3F20.8)") ts+1, tc+ dtini_g/60.0, i, j, q_g(i,j), qpx_g(i,j), qlatf                                     dx_ar_g(i,j)
!                        enddo


                else !* if (nldxj==ncomp-1) then
                !* Modify the current lengths of the segments of reach j so as to meet being larger than delxj_mn
                  !* modify_dx-step1: compute a new uniform segment length that are larger than delxj_mn -> delxj_m
                    sumdxj=0.0
                    do i=1, ncomp-1
                        sumdxj= sumdxj + dx_ar_g(i,j)
                    end do

                    mj= sumdxj/delxj_mn
                    ncomp_m= mj+1

                    if (ncomp_m>ncomp) then
                        dmyr= ncomp-1
                        delxj_m= sumdxj/dmyr
                        ncomp_m= ncomp
                    else
                        if (mj==0) then
                            delxj_m= sumdxj
                            ncomp_m=2
                        else
                            delxj_m= sumdxj/real(mj,KIND(sumdxj))
                            ncomp_m= ncomp_m
                        end if
                    end if

!                             write(13,"(3A10, 3A20)") "ts", "i", "j", "delxj_mn", "delxj_m",  "dx_ar_g(i,j)"
!                            do i=1, ncomp-1
!                                write(13,"(3I10, 2F20.2, F20.2)") ts, i, j, delxj_mn, delxj_m, dx_ar_g(i,j)
!                            enddo

                  !* modify_dx-step2: compute dx_ar at time ts for the changed length delxj_m
                    allocate (dx_ar_m_g(ncomp_m))
                    do i_m=1, ncomp_m-1
                        dx_ar_m_g(i_m)= delxj_m
                    enddo
                    dx_ar_m_g(ncomp_m)= dx_ar_m_g(ncomp_m-1) !* nominal value, not used.

                  !* modify_dx-step3: compute celty, diffty, q, and qpx at time ts for the changed length delxj_m
                    allocate (celty_m_g(ncomp_m), diffty_m_g(ncomp_m), q_m_g(ncomp_m), qpx_m_g(ncomp_m))
                    call mod_cdqqpx(j, ncomp, ncomp_m, delxj_m)

!                            write(13,"(3A10, 4A20)") "ts", "i", "j", "celty_g(i,j)", "diffty_g(i,j)", "q_g(i,j)", "qpx_g(i,j)"
!                            do i=1, ncomp
!                                write(13,"(3I10, 4F20.4)") ts, i, j, celty_g(i,j), diffty_g(i,j), q_g(i,j), qpx_g(i,j)
!                            enddo
!                            write(13,"(3A10, 4A20)") "ts", "i_m", "j", "celty_m_g(i_m)", "diffty_m_g(i_m)",&
!                                                        "q_m_g(i_m)", "qpx_m_g(i_m)"
!                            do i=1, ncomp_m
!                                write(13,"(3I10, 4F20.4)") ts, i, j, celty_m_g(i), diffty_m_g(i), q_m_g(i), qpx_m_g(i)
!                            enddo

                  !* modify_dx-step4: compute qlatj at time ts for the changed length delxj_m
                    allocate(qlatj_m_g(ncomp_m))
                    !* initialize
                    qlatj_m_g= initval
                    call mod_qlatj(j, ncomp, ncomp_m, delxj_m)
                            !* mass conservation test
                            dmyr1=0.0
                            dmyr2=0.0
                            do i=1, ncomp-1
                                dmyr1= dmyr1 + qlatj_g(i)*dx_ar_g(i,j)
                                dmyr2= dmyr2 + dx_ar_g(i,j)
                            end do
                            dmyr3=0.0
                            dmyr4=0.0
                            do i=1, ncomp_m-1
                                dmyr3= dmyr3 + qlatj_m_g(i)*dx_ar_m_g(i)
                                dmyr4= dmyr4 + dx_ar_m_g(i)
                            end do

!                            write(13,"(3A10, 4A20)") "ts", "i", "j", "qlatj_g(i)", "dx_ar_g(i,j)", "total_dist", "total_mass"
!                            do i=1, ncomp
!                                write(13,"(3I10, F20.8, 3F20.8)") ts, i, j, qlatj_g(i), dx_ar_g(i,j), dmyr2, dmyr1
!                            enddo
!                            write(13,"(3A10, 5A20)") "ts", "i_m", "j", "qlatj_m_g(i_m)", "dx_ar_m_g(i_m)", "total_dist",&
!                                                        "total_mass", "del_totalmass"
!                            do i=1, ncomp_m
!                                write(13,"(3I10, F20.8, 4F20.8)") ts, i, j, qlatj_m_g(i), dx_ar_m_g(i) ,dmyr4, dmyr3, dmyr1-dmyr3
!                            enddo

                  !* modify_dx-step5: Compute eei, ffi, exi, fxi for changed length delxj_m
                    allocate(eei_m_g(ncomp_m), ffi_m_g(ncomp_m), exi_m_g(ncomp_m), fxi_m_g(ncomp_m))
                    eei_m_g(1) = 1.0
                    ffi_m_g(1) = 0.
                    exi_m_g(1) = 0.
                    fxi_m_g(1) = 0.

                    call ef_calc_m(ncomp_m)
!
!                            write(91,"(4A10, 10A20)") "ts","CURR:tc_m", "i_m", "j", "eei_m_g", "ffi_m_g","exi_m_g","fxi_m_g",&
!                                                "q_m_g", "qpx_m_g", "qlatj_m_g", "celty_m_g", "diffty_m_g", "dx_ar_m_g"
!                            do i_m=1,ncomp_m
!                                write(91,"(I10, F10.1, 2I10, 10F20.8)") ts, tc, i_m, j,&
!                                                eei_m_g(i_m), ffi_m_g(i_m), exi_m_g(i_m), fxi_m_g(i_m),&
!                                                q_m_g(i_m), qpx_m_g(i_m), qlatj_m_g(i_m),&
!                                                celty_m_g(i_m), diffty_m_g(i_m), dx_ar_m_g(i_m)
!                            enddo

                   !* modify_dx-step6: Compute q and qpx at time ts+1 for changed length delxj_m
                    !* first, load lateral flow data at ts+1 instead of ts into qlatj_g
                    do i=1,ncomp-1
                        do n=1,nts_ql_g
                            !tarr_ql(n)= dt_ql_g*real(n,KIND(dt_ql_g))/60.0 !* [min]
                            varr_ql(n)= qlat_g(n,i,j) !* qlat_g(n,i,j) in unit of m2/sec
                        enddo
                        tf0= tc+dtini_g/60.0
                        qlatj_g(i)= intp_y(nts_ql_g, tarr_ql, varr_ql, tf0)
                    enddo

                    call mod_qlatj(j, ncomp, ncomp_m, delxj_m)

                    qlatf_m= qlatj_m_g(ncomp_m-1)
                    !* second, update q_m_g and qpx_m_g at ts+1
                    q_m_g(ncomp_m)=q_m_g(ncomp_m-1) + qlatf_m*dx_ar_m_g(ncomp_m-1) !* qlatf at ts+1 in unit of m2/sec
                    qpx_m_g(ncomp_m)=0.0

                    do i_m= ncomp_m-1, 1, -1
                        q_m_g(i_m)= eei_m_g(i_m) * q_m_g(i_m+1) + ffi_m_g(i_m)    !* q_g(i) at time t+dtini/60 [min] (=n+1)
                        qpx_m_g(i_m)= exi_m_g(i_m) *qpx_m_g(i_m+1) + fxi_m_g(i_m)  !* Qx(i) at time t+dtini/60 [min] (=n+1)
                    enddo
!
!                            write(91,"(4A10, 3A20)") "ts+1","FUR:tc_m", "i_m", "j", "q_m_g", "qpx_m_g", "qlatf_m"
!                            do i_m=1,ncomp_m
!                                write(91,"(I10, F10.1, 2I10, 3F20.8)") ts+1, tc+dtini_g/60.0, i_m, j,&
!                                                q_m_g(i_m), qpx_m_g(i_m), qlatf_m
!                            enddo

                  !* modify_dx-step7: Map q_m_g and qpx_m_g at ts+1 back to q_g and qpx_g at ts+1.  q_g at tst+1 becomes
                  !* input to computed elv, celty, diffty at ts+1.  qpx_g at ts+1 is also used in the next run of ef_calc
                  !* or ef_calc_m.
                    call mapback_qqpx(j, ncomp, ncomp_m, delxj_m)
!
!                            write(13,"(3A10, 2A20)") "ts", "i_m", "j", "q_m_g(i_m)_ts+1", "qpx_m_g(i_m)_f"
!                            do i=1, ncomp_m
!                                write(13,"(3I10, 2F20.4)") ts, i, j, q_m_g(i), qpx_m_g(i)
!                            enddo
!                            write(13,"(3A10, 2A20)") "ts", "i", "j", "q_g(i,j)_ts+1", "qpx_g(i,j)_ts+1"
!                            do i=1, ncomp
!                                write(13,"(3I10, 2F20.4)") ts, i, j, q_g(i,j), qpx_g(i,j)
!                            enddo
!
!                            write(91,"(4A10, 2A20)") "ts+1","FUR:tc_m", "MAPBACK:i", "j", "q_g", "qpx_g"
!                            do i=1,ncomp
!                                write(91,"(I10, F10.1, 2I10, 2F20.8)") ts+1, tc+dtini_g/60.0, i, j,&
!                                                q_g(i,j), qpx_g(i,j)
!                            enddo


                    deallocate (dx_ar_m_g)
                    deallocate (celty_m_g, diffty_m_g, q_m_g, qpx_m_g)
                    deallocate (qlatj_m_g)
                    deallocate(eei_m_g, ffi_m_g, exi_m_g, fxi_m_g)
                end if !* if (nldxj==ncomp-1) then ... else ...

                !* ------------------------------------------
                !* upstream boundary condition for  reach j.
                !* ------------------------------------------
                if (frnw_g(j,3)==0) then !* frnw_g(j,3) indicates the number of upstream reaches.
                !* head water reach
                    do n=1,nts_ub_g
                        !tarr_ub(n)= dt_ub_g*real(n,KIND(dt_ub_g))/60.0 !* [min]
                        varr_ub(n)= ubcd_g(n,j) !qbd(n,1,j)
                    enddo
                    tf0= tc+dtini_g/60.0
                    q_g(1,j)= intp_y(nts_ub_g, tarr_ub, varr_ub, tf0)
                    !** No zero discharge allowed
                    if (q_g(1,j)<0.0001) then !*<- 0.0001 is arbitrary value
                        q_g(1,j)= q_g(2,j)
                    end if
                else
                !* update q of ts+1 at a junction
                    qjt= 0.0
                    nusrch= frnw_g(j,3)  !* then number of upstream reaches
                    do rch=1, nusrch
                        usrchj= frnw_g(j,3+rch)  !* js corresponding to upstream reaches
                        ncomp_usrchj= frnw_g(usrchj,1)
                        qjt= qjt + q_g(ncomp_usrchj, usrchj)
                    enddo
                    q_g(1,j)= qjt
                endif

!                    write(91,"(4A10, A20)") "ts+1 ", "FUR:tc_M", "upbd:i=1", "j", "q_g"
!                    i=1
!                    write(91,"(I10,F10.1,2I10, F20.8)") ts+1, tc+dtini_g/60.0, i, j, q_g(i,j)
                deallocate(qlatj_g)

            enddo !* do j=1, nrch_g
            !write(91,*)

            !*---------------------------------------------------------------------------------------
            !* Compute y(=elv), celerity, diffusivity at ts+1 from downstream reach to upstream reach while
            !* computing from downstream node to upstream node for each reach.
            !*---------------------------------------------------------------------------------------
            do j=nrch_g, 1, -1

                 ncomp= frnw_g(j,1)  !ncomp= nx1(j)
                !* downstream boundary condition for elv at ts+1
                if (frnw_g(j,2)<0.0) then
                !* downstream boundary node
                    do n=1,nts_db_g
                        !tarr_db(n)= dt_db_g*real(n,KIND(dt_db_g))/60.0 !* [min]
                        varr_db(n)= dbcd_g(n)
                    enddo
                    tf0= tc+dtini_g/60.0
                    elv_g(ncomp,j)= intp_y(nts_db_g, tarr_db, varr_db, tf0)
                else
                !* downstream end node of a reach at a junction
                    dsrchj= frnw_g(j,2)    !* reach j's downstream reach index

                    elv_g(ncomp,j)= elv_g(1,dsrchj)
                endif

                !** Compute elevation, celerity, and diffusivity.
                call elv_calc(ts, ncomp, j)

!                write(91,"(4A10, 4A20)") "ts", "  tc", "i", "j","elv(ts+1,i,j)", "cel(ts+1,i,j)", "diff(ts+1,i,j)",&
!                                        "z_ar_g(i,j)"
!                do i=1, ncomp
!                    write(91,"(I10,F10.1,2I10, 4F20.8)") ts, tc, i, j, elv_g(i,j), celty_g(i,j),diffty_g(i,j),&
!                                                        z_ar_g(i,j)
!                enddo
            enddo !* j=nrch_g, 1, -1

            !cel_mx=maxval(cel_av_g)
            tc= tc + dtini_g/60.0 !* [min]
            ts= ts+1
            !* saveInterval in [sec]
            saveInterval_min=saveInterval_ev_g/60.0 !* [min]
            rmnd= mod(tc,saveInterval_min)  !* computes the remainder of the division of tc by saveinterval_min

            if ( (rmnd<0.0001).or.(tc==tfin_g*60.0) ) then
                do j=1, nrch_g
                    ncomp= frnw_g(j,1)
                    do i=1, ncomp
                        q_ev_g(ts_ev, i, j)=q_g(i,j)
                        elv_ev_g(ts_ev, i, j)= elv_g(i,j)
                        write(1,"(f10.3, 3I10, 3F20.4)") tc, ts_ev, i, j, q_ev_g(ts_ev,i,j), elv_ev_g(ts_ev,i,j),&
                                                            elv_ev_g(ts_ev,i,j)-z_ar_g(i,j)
                        !write(*,"(I5,F8.1,I2,I4, 2F20.2)") ts, tc,i,j, q_g(i,j), elv_g(i,j)-z_ar_g(i,j)
                    enddo
                enddo
                ts_ev=ts_ev+1
            endif

        enddo !*do while ( tc .lt. tfin *60.0)

        !deallocate(celty_g, diffty_g, q, elv_g, qpx_g, q_ev_g, elv_ev_g)
        !deallocate(tarr_ql, varr_ql, tarr_ub, varr_ub, tarr_db, varr_db, cel_av_g)
        !deallocate(ufhlt_m_g, ufqlt_m_g, ufhlt_f_g, ufqlt_f_g)
        !deallocate(frnw_g)
        !deallocate(qlat_g, ubcd_g, dbcd_g)
        !deallocate(z_ar_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g, So_ar_g)
        !deallocate(mann_ar_g, manncc_ar_g, dx_ar_g)

        deallocate(celty_g, diffty_g)
        deallocate(q_g, elv_g, qpx_g)
        !deallocate(qlat_g)
        !deallocate(q_ev_g, elv_ev_g)
        deallocate(tarr_ql, varr_ql)
        deallocate(tarr_ub, varr_ub)
        deallocate(tarr_db, varr_db)
        deallocate(cel_av_g)
        !deallocate(ufhlt_m_g, ufqlt_m_g)
        !deallocate(ufhlt_f_g, ufqlt_f_g)
        !deallocate(frnw_g)
        !deallocate( z_ar_g, bo_ar_g, traps_ar_g)
        !deallocate( tw_ar_g, twcc_ar_g )
        !deallocate( mann_ar_g, manncc_ar_g )
        !deallocate( So_ar_g, dx_ar_g )
        !deallocate( dbcd_g )
        !deallocate( ubcd_g )

    endsubroutine routenw

!nd program main2
endmodule mdiff

