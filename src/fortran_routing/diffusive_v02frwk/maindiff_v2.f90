module mdiffv2

    implicit none

    double precision, dimension(:,:), allocatable :: celty_g, diffty_g, q_g, elv_g, qpx_g
    double precision, dimension(:), allocatable :: cel_av_g
    double precision, dimension(:), allocatable :: qlatj_g
    double precision, dimension(:), allocatable :: eei_g, ffi_g, exi_g, fxi_g
    double precision, dimension(:), allocatable :: celty_m_g, diffty_m_g, q_m_g, qpx_m_g
    double precision, dimension(:), allocatable :: qlatj_m_g
    double precision, dimension(:), allocatable :: eei_m_g, ffi_m_g, exi_m_g, fxi_m_g
    double precision, dimension(:), allocatable :: dx_ar_m_g
    double precision, dimension(:,:,:), allocatable :: bfuf_g    
    double precision:: z_g, bo_g, traps_g, tw_g, twcc_g, So_g, mann_g, manncc_g


contains
    !*--------------------------------------------------------------------------------
    !*       Route network by using Python-to-Fortran network traversal map together
    !*          with diffusive routing engines
    !
    !*--------------------------------------------------------------------------------
    subroutine diffnw(dtini_g, t0_g, tfin_g, saveinterval_g, saveinterval_ev_g, dt_ql_g, dt_ub_g, dt_db_g, &
                        nts_ql_g, nts_ub_g, nts_db_g, &
                        mxncomp_g, nrch_g, z_ar_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g, &
                        mann_ar_g, manncc_ar_g, so_ar_g, dx_ar_g, &
                        nhincr_m_g, nhincr_f_g, ufhlt_m_g,  ufqlt_m_g, ufhlt_f_g, ufqlt_f_g, &
                        frnw_col, dfrnw_g, qlat_g, ubcd_g, dbcd_g, &
                        cfl_g, theta_g, tzeq_flag_g, y_opt_g, so_llm_g, &
                        ntss_ev_g, q_ev_g, elv_ev_g)

        implicit none

        integer, intent(in) :: mxncomp_g, nrch_g
        integer, intent(in) :: nts_ql_g, nts_ub_g, nts_db_g, ntss_ev_g
        integer, intent(in) :: nhincr_m_g, nhincr_f_g, frnw_col
        double precision,intent(in) :: dtini_g, t0_g, tfin_g, saveinterval_g, saveinterval_ev_g, dt_ql_g, dt_ub_g, dt_db_g
        double precision, dimension(mxncomp_g, nrch_g), intent(in) :: z_ar_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g
        double precision, dimension(mxncomp_g, nrch_g), intent(in) :: mann_ar_g, manncc_ar_g, dx_ar_g
        double precision, dimension(mxncomp_g, nrch_g, nhincr_m_g), intent(in) :: ufhlt_m_g,  ufqlt_m_g
        double precision, dimension(mxncomp_g, nrch_g, nhincr_f_g), intent(in) :: ufhlt_f_g, ufqlt_f_g
        double precision, dimension(nrch_g, frnw_col), intent(in) :: dfrnw_g !* set frnw_col=8
        double precision, dimension(nts_ql_g, mxncomp_g, nrch_g), intent(in) :: qlat_g
        double precision, dimension(nts_ub_g, nrch_g), intent(in) :: ubcd_g
        double precision, dimension(nts_db_g), intent(in) :: dbcd_g 

        double precision, intent(in) :: cfl_g, theta_g, so_llm_g
        integer, intent(in) :: tzeq_flag_g !* 0 for lookup tabale; 1 for using procedures to compute trapz.ch.geo.
        integer, intent(in) :: y_opt_g  !* 1 for normal depth(kinematic); 2 for dept of diffusive wave.
        double precision, dimension(mxncomp_g, nrch_g), intent(inout) :: so_ar_g
        double precision, dimension(ntss_ev_g, mxncomp_g, nrch_g), intent(out) :: q_ev_g, elv_ev_g
        ! Local storage
        integer :: ts, n, i, j, ts_ev, i1
        integer :: ncomp
        double precision ::  tc, tf0
        double precision :: qlatf
        integer :: nusrch, rch, usrchj, ncomp_usrchj, dsrchj
        integer :: lntss, lmxncomp, lnrch
        double precision :: qjt
        double precision, dimension(:), allocatable ::tarr_ql, varr_ql, tarr_ub, varr_ub, tarr_db, varr_db
        integer :: mi
        double precision :: vr, vr2, mtp, mtp2, mn, mn2, yintc, cel_mx, rmnd, saveInterval_min
        double precision :: dmyr, dmyr1, dmyr2, dmyr3, dmyr4 !test only
        integer :: nldxj, mj, ncomp_m, i_m, x_mflag, i_s, i2, xflag, i_m_s
        double precision ::  delxj_mn, delxj_m, sumdxj, x_m, sumx1, sumx2, u1, u2, u_m
        double precision :: sumqlatdx, sumqlatx_m, sumqlatdx_m, qlatf_m, x_o, sumx1_m, sumx2_m, u_o
        double precision :: so_us, c_us, so_ds, c_ds, adjso, initval
        integer :: frj, iseg
        integer, dimension(:,:), allocatable :: frnw_g
        allocate(frnw_g(nrch_g, frnw_col))
        frnw_g=int(dfrnw_g) 

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
!*Py         nts_ql_g= int( (tfin_g - t0_g)*3600.0/dt_ql_g+1, KIND(nts_ql_g) ) !* based on time of [sec]!*Py         nts_ub_g= nts_ql_g
!*Py         nts_db_g= int( (tfin_g - t0_g)*3600.0/dt_db_g+1, KIND(nts_db_g) ) !* based on time of [sec]
        allocate(celty_g(mxncomp_g, nrch_g), diffty_g(mxncomp_g, nrch_g))
        allocate(q_g(mxncomp_g, nrch_g), elv_g(mxncomp_g, nrch_g), qpx_g(mxncomp_g, nrch_g))
        allocate(tarr_ql(nts_ql_g), varr_ql(nts_ql_g))
        allocate(tarr_ub(nts_ub_g), varr_ub(nts_ub_g))
        allocate(tarr_db(nts_db_g), varr_db(nts_db_g))
        allocate(cel_av_g(nrch_g))

        !++-----------------------------------------------------------------------------------+
        !+        adjust abnormally small channel bottom slope
        !+ Definition of abnormal slope:
        !+  1. slopes less than a chosen lower limit (=so_llm)
        !+ Repair method:
        !+  1. take average of slopes of adjacent segments
        !++-----------------------------------------------------------------------------------+
!*Py    so_llm_g=0.0005 !*initial guess
        do j=1,nrch_g
            ncomp=frnw_g(j,1)
            do i=1, ncomp-1
                if (so_ar_g(i,j).lt.so_llm_g) then
                    !* adjacent upstream segment's slope
                    so_us=0.0
                    c_us=0.0
                    i2=i-1
                    do while (i2.ge.1)
                        if (so_ar_g(i2,j).ge.so_llm_g) then
                            so_us= so_ar_g(i2,j)
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
                        if (so_ar_g(i2,j).ge.so_llm_g) then
                            so_ds= so_ar_g(i2,j)
                            c_ds=1.0
                            exit
                        endif
                        i2=i2+1
                    end do
                    if (c_us+c_ds.gt.0.0) then
                        adjso= (so_us+so_ds)/(c_us+c_ds)
                        so_ar_g(i,j)= adjso
                    else
                        so_ar_g(i,j)= so_llm_g
                    endif
                endif
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
        allocate( bfuf_g(mxncomp_g, nrch_g,2) )
        call bfdep_ufQ(mxncomp_g, nrch_g, z_ar_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g,&
                        mann_ar_g, manncc_ar_g, dx_ar_g, so_ar_g, frnw_col, frnw_g)


        !++-------------------------------------------+
        !+              Lateral inflow
        !+
        !++-------------------------------------------
!*Py      do j=1, nrch_g
!            ncomp=frnw_g(j,1)
!            do i=1, ncomp
!                do n=1, nts_ql_g
!                    read(31,*) frj, iseg, tarr_ql(n), qlat_g(n,i,j) !* qlat_g(n,i,j) in [m^2/sec]
!                end do
!            end do
!        enddod
        !* time step series for lateral flow
        do n=1, nts_ql_g
            tarr_ql(n)= t0_g*60.0 + dt_ql_g*real(n-1,KIND(dt_ql_g))/60.0 !* [min]
        end do

        !++--------------------------------------------------------------+
        !+              Upstream boundary condition
        !+
        !+ Assume that hydrographs are always available at the upper ends
        !+ of the most upstream links.
        !++--------------------------------------------------------------+
!*Py     allocate( ubcd_g(nts_ub_g, nrch_g) )
!        do j=1, nrch_g
!            do n=1,nts_ub_g
!                ubcd_g(n,j)=-1.0 !* [m^3/sec]
!            end do
!        end do
        !* time step series for upstream boundary data
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
!        do n=1, nts_db_g
!            read(33,*) tarr_db(n), dbcd_g(n) !*  dbcd_g(n)  in [meter] based on NAD88 meter datum that RouteLink uses
!        end do
        !* time step series for downstream boundary data
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

        tc = t0_g*60.0     !!! t0 is in hour. tc is in minutes
        ts=1    !*simulation time step in parallel with tc.
        ts_ev=1 !* time step for outputting q and elv at evaluation time
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
                delxj_mn= cel_av_g(j)*dtini_g/cfl_g
                ncomp= frnw_g(j,1)
                !* interpolate qlat_g at the current time tc
                allocate(qlatj_g(ncomp))
                qlatj_g=initval !* initialize elements' dummy values
                do i=1,ncomp-1
                    do n=1,nts_ql_g
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
                    call ef_calc(mxncomp_g, nrch_g, ncomp, j, dtini_g, theta_g, &
                                dx_ar_g)
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
                  !* modify_dx-step2: compute dx_ar at time ts for the changed length delxj_m
                    allocate (dx_ar_m_g(ncomp_m))
                    do i_m=1, ncomp_m-1
                        dx_ar_m_g(i_m)= delxj_m
                    enddo
                    dx_ar_m_g(ncomp_m)= dx_ar_m_g(ncomp_m-1) !* nominal value, not used.
                  !* modify_dx-step3: compute celty, diffty, q, and qpx at time ts for the changed length delxj_m
                    allocate (celty_m_g(ncomp_m), diffty_m_g(ncomp_m), q_m_g(ncomp_m), qpx_m_g(ncomp_m))
                    !call mod_cdqqpx(j, ncomp, ncomp_m, delxj_m)
                    call mod_cdqqpx(mxncomp_g, nrch_g, j, ncomp, ncomp_m, dx_ar_g)
                  !* modify_dx-step4: compute qlatj at time ts for the changed length delxj_m
                    allocate(qlatj_m_g(ncomp_m))
                    !* initialize
                    qlatj_m_g= initval
                    !call mod_qlatj(j, ncomp, ncomp_m, delxj_m)
                    call mod_qlatj(mxncomp_g, nrch_g, j, ncomp, ncomp_m, dx_ar_g, delxj_m)
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
                  !* modify_dx-step5: Compute eei, ffi, exi, fxi for changed length delxj_m
                    allocate(eei_m_g(ncomp_m), ffi_m_g(ncomp_m), exi_m_g(ncomp_m), fxi_m_g(ncomp_m))
                    eei_m_g(1) = 1.0
                    ffi_m_g(1) = 0.
                    exi_m_g(1) = 0.
                    fxi_m_g(1) = 0.
                    call ef_calc_m(ncomp_m, dtini_g, theta_g)
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

                    call mod_qlatj(mxncomp_g, nrch_g, j, ncomp, ncomp_m, dx_ar_g, delxj_m)

                    qlatf_m= qlatj_m_g(ncomp_m-1)
                    !* second, update q_m_g and qpx_m_g at ts+1
                    q_m_g(ncomp_m)=q_m_g(ncomp_m-1) + qlatf_m*dx_ar_m_g(ncomp_m-1) !* qlatf at ts+1 in unit of m2/sec
                    qpx_m_g(ncomp_m)=0.0

                    do i_m= ncomp_m-1, 1, -1
                        q_m_g(i_m)= eei_m_g(i_m) * q_m_g(i_m+1) + ffi_m_g(i_m)    !* q_g(i) at time t+dtini/60 [min] (=n+1)
                        qpx_m_g(i_m)= exi_m_g(i_m) *qpx_m_g(i_m+1) + fxi_m_g(i_m)  !* Qx(i) at time t+dtini/60 [min] (=n+1)
                    enddo
                  !* modify_dx-step7: Map q_m_g and qpx_m_g at ts+1 back to q_g and qpx_g at ts+1.  q_g at tst+1 becomes
                  !* input to computed elv, celty, diffty at ts+1.  qpx_g at ts+1 is also used in the next run of ef_calc
                  !* or ef_calc_m.
                    call mapback_qqpx(mxncomp_g, nrch_g, j, ncomp, ncomp_m, dx_ar_g, delxj_m)

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
                        varr_ub(n)= ubcd_g(n,j)
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
                deallocate(qlatj_g)
            enddo !* do j=1, nrch_g
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
                call elv_calc(ts, ncomp, j, tzeq_flag_g, y_opt_g, &
                        mxncomp_g, nrch_g, z_ar_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g,&
                        mann_ar_g, manncc_ar_g, dx_ar_g, &
                        nhincr_m_g, nhincr_f_g, ufhlt_m_g,  ufqlt_m_g, ufhlt_f_g, ufqlt_f_g)
            enddo !* j=nrch_g, 1, -1

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
                    enddo
                enddo
                ts_ev=ts_ev+1
            endif

        enddo !*do while ( tc .lt. tfin *60.0)

        deallocate(celty_g, diffty_g)
        deallocate(q_g, elv_g, qpx_g)
        deallocate(tarr_ql, varr_ql)
        deallocate(tarr_ub, varr_ub)
        deallocate(tarr_db, varr_db)
        deallocate(cel_av_g)
	deallocate(frnw_g)

    endsubroutine diffnw

!*--------------------------------------------------------------------------------
!*          Compute eei_g, ffi_g, exi_g, and fxi_g from upstream to downstream
!
!       input: dtini, dx_ar_g, celty_g, diffusivity, q, qpx_g, qlat_g
!
! %%%%% qlat_g takes a unit of m^2/sec
!*--------------------------------------------------------------------------------
    subroutine ef_calc(mxncomp_g, nrch_g, ncomp, j, dtini_g, theta_g, dx_ar_g)
        implicit none

        integer, intent(in) :: mxncomp_g, nrch_g, ncomp, j
        double precision, intent(in) :: dtini_g, theta_g
        double precision, dimension(mxncomp_g, nrch_g), intent(in) :: dx_ar_g

        integer :: i
        double precision :: cour !,  theta
        double precision :: a1, a2, a3, a4, b1, b2, b3, b4, dd1, dd2, dd3, dd4, h1, h2, h3, h4
        double precision :: qy, qxy, qxxy, qxxxy, ppi, qqi, rri, ssi, sxi, cour2

        do i = 2,ncomp
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

            eei_g(i) = -1.0 * rri / ( ppi * eei_g(i-1) + qqi )
            ffi_g(i) = ( ssi - ppi * ffi_g(i-1) ) / ( ppi * eei_g(i-1) + qqi )
            exi_g(i) = -1.0 * rri / ( ppi * exi_g(i-1) + qqi )
            fxi_g(i) = ( sxi - ppi * fxi_g(i-1) ) / ( ppi * exi_g(i-1) + qqi )
        enddo

    end subroutine ef_calc
!*-------------------------------------------------------------------------------------
!*          Compute eei_m_g, ffi_m_g, exi_m_g, and fxi_m_g from upstream to downstream
!*          for modified segment lengths of a given reach j
!
!       input: dtini, dx_ar_m_g, celty_m_g, diffty_m_g, q_m_g, qpx_m_g, qlat_m_g
!
! %%%%% qlat_m_g takes a unit of m^2/sec
!*--------------------------------------------------------------------------------------
    subroutine ef_calc_m(ncomp_m, dtini_g, theta_g)

        implicit none

        integer, intent(in) :: ncomp_m
        double precision, intent(in) :: dtini_g, theta_g

        integer :: i
        double precision :: cour
        double precision :: a1, a2, a3, a4, b1, b2, b3, b4, dd1, dd2, dd3, dd4, h1, h2, h3, h4
        double precision :: qy, qxy, qxxy, qxxxy, ppi, qqi, rri, ssi, sxi, cour2

        do i = 2, ncomp_m
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

            ssi = qy  + dtini_g * diffty_m_g(i) * ( 1.0 - theta_g ) * qxxy + dtini_g * celty_m_g(i) * qlatj_m_g(i-1)
            sxi = qxy + dtini_g * diffty_m_g(i) * ( 1.0 - theta_g ) * qxxxy+ dtini_g * celty_m_g(i) * qlatj_m_g(i-1)/ dx_ar_m_g(i-1)

            eei_m_g(i) = -1.0 * rri / ( ppi * eei_m_g(i-1) + qqi )
            ffi_m_g(i) = ( ssi - ppi * ffi_m_g(i-1) ) / ( ppi * eei_m_g(i-1) + qqi )
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
    subroutine elv_calc(ts, ncomp, j, tzeq_flag_g, y_opt_g, &
                        mxncomp_g, nrch_g, z_ar_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g,&
                        mann_ar_g, manncc_ar_g, dx_ar_g, &
                        nhincr_m_g, nhincr_f_g, ufhlt_m_g,  ufqlt_m_g, ufhlt_f_g, ufqlt_f_g)

        implicit none

        integer, intent(in) :: mxncomp_g, nrch_g, ts, ncomp, j
        integer, intent(in) :: tzeq_flag_g, y_opt_g
        integer, intent(in) :: nhincr_m_g, nhincr_f_g
        double precision, dimension(mxncomp_g, nrch_g), intent(in) :: z_ar_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g
        double precision, dimension(mxncomp_g, nrch_g), intent(in) :: mann_ar_g, manncc_ar_g, dx_ar_g
        double precision, dimension(mxncomp_g, nrch_g, nhincr_m_g), intent(in) :: ufhlt_m_g,  ufqlt_m_g
        double precision, dimension(mxncomp_g, nrch_g, nhincr_f_g), intent(in) :: ufhlt_f_g, ufqlt_f_g

        integer :: i, icol, i1
        double precision ::  q_sk_multi, sfi, xt , width
        double precision :: depth_tz, area_tz, hydr_tz, conv_tz, hbf_tz, ybf_tz, topwd_tz, emann_tz, sum1, sum2, celav, diffav
        double precision, allocatable ::  co(:)
        double precision :: dsc, hnm
        double precision, dimension(:), allocatable ::  harr_m, qarr_m, harr_f, qarr_f
        double precision ::  hbf, ufQbf
        allocate(co(ncomp))
        allocate(harr_m(nhincr_m_g), qarr_m(nhincr_m_g), harr_f(nhincr_f_g), qarr_f(nhincr_f_g))
        !*-------------------------------------------------------------------------------
        !*    Partial diffusive elevation(=normal depth)
        !*    elv at i <- normal_depth_LT{q at i, ch.geo at i}, where LT is lookup table
        !*-------------------------------------------------------------------------------
        if (y_opt_g==1) then
            do i= ncomp-1, 1, -1 !* by p.3,RM5
                dsc= abs(q_g(i,j))
                !* By lookup tables
                if (tzeq_flag_g==0) then
                !* lookup table created from x-section approximation.
!                    do i1=1, nel_g
!                        elvarr(i1)= xsec_attr_rch_g(i, j, i1, 1)   !* elevation  i1=1,...,nel
!                        qarr(i1)= xsec_attr_rch_g(i, j, i1, 4)     !* uniform flow
!                    enddo
!                    elv_g(i,j)= intp_y(nel_g, qarr, elvarr, dsc)
                elseif (tzeq_flag_g==1) then
                !* lookup table created from x-section equations for trap.main & rect.floodplains
                    !if (dsc<= ufqlt_m_g(i, j, nhincr_m_g)) then
                    !* inbank flow
                        !do i1=1,nhincr_m_g
                        !    qarr_m(i1)= ufqlt_m_g(i, j, i1)
                        !    harr_m(i1)= ufhlt_m_g(i, j, i1)
                        !enddo
                        !hnm= intp_y(nhincr_m_g, qarr_m, harr_m, dsc)
                     !else
                        !* overbank flow
                        !do i1=1,nhincr_f_g
                        !    qarr_f(i1)= ufqlt_f_g(i, j, i1)
                        !    harr_f(i1)= ufhlt_f_g(i, j, i1)
                        !enddo
                        !hnm= intp_y(nhincr_f_g, qarr_f, harr_f, dsc)
                    !endif
                    z_g= z_ar_g(i,j)
                    bo_g= bo_ar_g(i,j)
                    traps_g= traps_ar_g(i,j)
                    tw_g= tw_ar_g(i,j)
                    twcc_g= twcc_ar_g(i,j)
                    mann_g= mann_ar_g(i,j)  !1.0/sk(i,j)
                    manncc_g= manncc_ar_g(i,j) ! 1.0/skCC1(i,j)
                    hbf= bfuf_g(i,j,1)      !* bankfull water depth
                    ufQbf= bfuf_g(i,j,2)    !* bankfull uniform flow

                    call bsec_nmdep(dsc, hbf, ufQbf, hnm)
                    elv_g(i,j)= hnm + z_ar_g(i,j)
                endif
            enddo !* do i= ncomp-1, 1, -1
        endif !* if (y_opt_g==1) then

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
!                icol= 2 !<-5
!                co(i)= q_sk_multi*intp_elev(icol, i, j, xt)
!                icol= 3 !<-6
!                width= intp_elev(icol, i, j, xt)
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
            ybf_tz= hbf_tz + z_g  !* bankfull depth/elevation for trapezoidal main channel
            if (xt.le.ybf_tz) then
            !* use Manning's N for main channel
                celty_g(i,j)= 5.0 / 3.0 * sfi ** 0.3 * abs(q_g(i,j)) ** 0.4 / width ** 0.4 / ( mann_g*(1/q_sk_multi)) ** 0.6
                emann_tz= -100.0 !* test only
            else
            !* when overbankflow, use equivalent manning's N
                depth_tz= xt - z_g
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
        enddo   !*do i=ncomp,2,-1

        !* averaged celty_g thorough all segments to prevent a wavefront of a cell with too high speed which
        !* produces wave breaking and instability
        celav= sum1/real(ncomp-1,KIND(sum1))
        diffav= sum2/real(ncomp-1,KIND(sum2))
        do i=1,ncomp
            celty_g(i,j) = celav
            diffty_g(i,j)= diffav
            if (diffty_g(i,j)>10.0) diffty_g(i,j)=10.0
        end do
        cel_av_g(j)=celav

        deallocate(co)
        deallocate(harr_m, qarr_m, harr_f, qarr_f)

    end subroutine elv_calc

!*-----------------------------------------------------------*!






!**              Linear Interpolation Tools                 **!






!*-----------------------------------------------------------*!
    !*--------------------------------------------
    !           Interpolate any value
    !
    !*--------------------------------------------
    double precision function intp_y(nrow, xarr, yarr, x)
        implicit none
        integer, intent(in) :: nrow
        double precision, dimension(nrow), intent(in) :: xarr, yarr
        double precision, intent(in) :: x
        integer :: irow
        double precision :: x1, y1, x2, y2, y

        irow= locate(xarr, x)
        if (irow.eq.0) irow= 1
        if (irow.eq.nrow) irow= nrow-1
        x1= xarr(irow); y1= yarr(irow)
        x2= xarr(irow+1); y2= yarr(irow+1)
        y= LInterpol(x1,y1,x2,y2,x)
        intp_y = y

    end function intp_y
    !*----------------------------------------------------
    !           Interpolate channel property by elev
    ! ** used only for elv_calc
    !*----------------------------------------------------
!    double precision function intp_elev(icol, inode, j, x)
!        implicit none
!        integer(KIND=i4b), intent(in) :: icol, inode, j
!        real(KIND=dp), intent(in) :: x
!        real(KIND=dp), allocatable :: dmyv(:)
!        integer(KIND=i4b) ::  i1, irow, nrow
!        double precision :: x1, y1, x2, y2, y
!
!        nrow= nel_g
!        allocate(dmyv(nrow))
!        do i1=1, nrow
!            dmyv(i1)= xsec_attr_rch_g(inode, j, i1, 1 )  !* elevation vector
!        end do
!        irow= locate(dmyv, x)
!        if (irow.eq.0) irow= 1
!        if (irow.eq.nrow) irow= nrow-1
!        x1= xsec_attr_rch_g(inode, j, irow, 1 )
!        y1= xsec_attr_rch_g(inode, j, irow, icol)
!        x2= xsec_attr_rch_g(inode, j, irow+1, 1)
!        y2= xsec_attr_rch_g(inode, j, irow+1, icol)
!
!        y= LInterpol(x1,y1,x2,y2,x)
!        intp_elev= y
!        deallocate(dmyv)
!
!    end function intp_elev
    !*-----------------------------------------------------------------------------
    !               Locate function in f90, p.1045,NR f90
    !
    !   klo=max(min(locate(xa,x),n-1),1) In the Fortran 77 version of splint,
    !   there is in-line code to find the location in the table by bisection. Here
    !   we prefer an explicit call to locate, which performs the bisection. On
    !   some massively multiprocessor (MMP) machines, one might substitute a different,
    !   more parallel algorithm (see next note).
    !*-----------------------------------------------------------------------------
    integer function locate(xx,x)

        implicit none
        double precision, dimension(:), intent(in) :: xx
        double precision, intent(in) :: x
        !* Given an array xx(1:N), and given a value x, returns a value j such that x is between
        !* xx(j) and xx(j + 1). xx must be monotonic, either increasing or decreasing.
        !* j = 0 or j = N is returned to indicate that x is out of range.
        integer :: n,jl,jm,ju
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
    double precision function LInterpol(x1,y1,x2,y2,x)

        implicit none
        double precision, intent(in) :: x1, y1, x2, y2, x

        !* interpolate y for the given x
        LInterpol= (y2-y1)/(x2-x1)*(x-x1)+y1

    end function LInterpol
!*-----------------------------------------------------------*!






!**                 Reach Tuning Tools                      **!






!*-----------------------------------------------------------*!
    !*--------------------------------------------------------------------------
    !* Compute qlatj at a given time for the changed length delxj_m
    !* Refer to p.2-5, RM5.
    !*--------------------------------------------------------------------------
    subroutine mod_qlatj(mxncomp_g, nrch_g, j, ncomp, ncomp_m, dx_ar_g, delxj_m)

        implicit none
        integer, intent(in) :: mxncomp_g, nrch_g, j, ncomp, ncomp_m
        double precision, intent(in) :: delxj_m
        double precision, dimension(mxncomp_g, nrch_g), intent(in) :: dx_ar_g

        integer :: i, i_m, i_s, i2, x_mflag
        double precision :: x_m, sumx1, sumx2, sumqlatdx, dx_s, sumqlatx_m, sumqlatdx_m

        x_m=0.0
        do i_m=1, ncomp_m-1
            x_m= x_m + dx_ar_m_g(i_m)
            x_mflag=0
            i=1
            sumx1=0.0
            sumx2=0.0
            do while ((x_mflag.eq.0).and.(i.le.ncomp-1))
                sumx1= sumx1 + dx_ar_g(i,j)
                sumx2= sumx1 + dx_ar_g(i+1,j)
                if (x_m.le.dx_ar_g(1,j)) then
                    i_s=1
                    sumx2= dx_ar_g(1,j)
                    x_mflag=1
                elseif ((x_m.gt.sumx1).and.(x_m.le.sumx2)) then
                    i_s= i+1
                    x_mflag=1
                endif
                i=i+1
            enddo

            sumqlatdx=0.0
            do i2= 1, i_s
                sumqlatdx= sumqlatdx + qlatj_g(i2)*dx_ar_g(i2,j)
            enddo
            dx_s= sumx2 - x_m
            !* sum of qlat up to distance x_m from node1(=distance zero)
            sumqlatx_m= sumqlatdx - qlatj_g(i_s)*dx_s
            !* sum of qlat_m up to index i_m - 1
            sumqlatdx_m=0.0
            if (i_m>1) then
                do i2=1, i_m-1
                    sumqlatdx_m = sumqlatdx_m + qlatj_m_g(i2)*delxj_m
                enddo
            endif
            qlatj_m_g(i_m)= (sumqlatx_m - sumqlatdx_m)/delxj_m
        enddo !* do i_m=1, ncomp_m-1

    end subroutine mod_qlatj

    !* ----------------------------------------------------------------------------
    !* compute celty, diffty, q, and qpx at time ts for the changed length delxj_m
    !* p.2-5, RM5
    !* ----------------------------------------------------------------------------
    subroutine mod_cdqqpx(mxncomp_g, nrch_g, j, ncomp, ncomp_m, dx_ar_g)

        implicit none
        integer,intent(in) :: mxncomp_g, nrch_g, j, ncomp, ncomp_m
        double precision, dimension(mxncomp_g, nrch_g), intent(in) :: dx_ar_g

        integer :: i, i_m, i_s, ifc
        double precision :: x_m, x_mflag, sumx1, sumx2, dx_s, u1, u2, u_m

        !* first node
        celty_m_g(1) = celty_g(1,j)
        diffty_m_g(1) = diffty_g(1,j)
        q_m_g(1) = q_g(1,j)
        qpx_m_g(1)= qpx_g(1,j)
        !* last node
        celty_m_g(ncomp_m) = celty_g(ncomp,j)
        diffty_m_g(ncomp_m) = diffty_g(ncomp,j)
        q_m_g(ncomp_m) = q_g(ncomp,j)
        qpx_m_g(ncomp_m)= qpx_g(ncomp,j)
        !* in-between nodes
        x_m=0.0
        do i_m=2, ncomp_m-1
            x_m= x_m + dx_ar_m_g(i_m-1) !delxj_m*real(i_m-1,KIND(delxj_m))
            x_mflag=0
            i=1
            sumx1=0.0
            sumx2=0.0
            do while ((x_mflag.eq.0).and.(i.le.ncomp-1))
                sumx1= sumx1 + dx_ar_g(i,j)
                sumx2= sumx1 + dx_ar_g(i+1,j)
                if (x_m.le.dx_ar_g(1,j)) then
                    i_s=1
                    sumx1=0.0
                    sumx2= dx_ar_g(1,j)
                    x_mflag=1
                elseif ((x_m.gt.sumx1).and.(x_m.le.sumx2)) then
                    i_s= i+1
                    x_mflag=1
                endif
                i=i+1
            enddo
            dx_s= x_m- sumx1
            !dx_i_s= dx_ar_g(i_s,j)
            !delx= dx_i_s - dx_s
            do ifc=1,4
                if (ifc==1) then
                    u1= celty_g(i_s,j)
                    u2= celty_g(i_s+1,j)
                elseif (ifc==2) then
                    u1= diffty_g(i_s,j)
                    u2= diffty_g(i_s+1,j)
                elseif (ifc==3) then
                    u1= q_g(i_s,j)
                    u2= q_g(i_s+1,j)
                elseif (ifc==4) then
                    u1= qpx_g(i_s,j)
                    u2= qpx_g(i_s+1,j)
                endif
                !u_m= (u2-u1)*(delx-dx_i_s)/dx_i_s + u2
                u_m= (u2-u1)*dx_s/dx_ar_g(i_s,j) + u1

                if (ifc==1) then
                    celty_m_g(i_m) = u_m
                elseif (ifc==2) then
                    diffty_m_g(i_m) = u_m
                elseif (ifc==3) then
                    q_m_g(i_m) = u_m
                elseif (ifc==4) then
                    qpx_m_g(i_m) = u_m
                endif
            enddo
        enddo !* do i_m=2, ncomp_m-1

    end subroutine mod_cdqqpx

    !* ----------------------------------------------------------------------------------
    !* Map q_m_g and qpx_m_g at ts+1 back to q_g and qpx_g at ts+1.  q_g at tst+1 becomes
    !* input to computed elv, celty, diffty at ts+1.  qpx_g at ts+1 is also used in the
    !* next run of ef_calc or ef_calc_m. p2-10,RM5
    !* ----------------------------------------------------------------------------
    subroutine mapback_qqpx(mxncomp_g, nrch_g, j, ncomp, ncomp_m, dx_ar_g, delxj_m)

        implicit none
        integer, intent(in) :: mxncomp_g, nrch_g, j, ncomp, ncomp_m
        double precision,intent(in) :: delxj_m
        double precision, dimension(mxncomp_g, nrch_g), intent(in) :: dx_ar_g

        integer :: xflag, i, i_m, i_m_s, ifc
        double precision :: x_o, sumx1_m, sumx2_m, dx_s, u1, u2, u_o  !delx_m,

        q_g(1,j)= q_m_g(1)
        qpx_g(1,j)= qpx_m_g(1)
        q_g(ncomp,j)= q_m_g(ncomp_m)
        qpx_g(ncomp,j)= qpx_m_g(ncomp_m)

        x_o= 0
        do i=2, ncomp-1
            x_o= x_o + dx_ar_g(i-1,j) !* x distance from upstream end node of original segment lengths.
            xflag=0
            i_m=0
            sumx1_m=0.0
            sumx2_m=0.0
            do while ((xflag.eq.0).and.(i_m.lt.ncomp_m))
                sumx1_m= delxj_m*real(i_m,KIND(delxj_m))
                sumx2_m= delxj_m*real(i_m+1,KIND(delxj_m))
                if ((x_o.gt.sumx1_m).and.(x_o.le.sumx2_m)) then
                    i_m_s= i_m+1
                    xflag=1
                endif
                i_m= i_m + 1
            enddo
            dx_s=  x_o - sumx1_m
            !delx_m= delxj_m - dx_s
            do ifc=1,2
                if (ifc==1) then
                    u1= q_m_g(i_m_s)
                    u2= q_m_g(i_m_s+1)
                elseif (ifc==2) then
                    u1= qpx_m_g(i_m_s)
                    u2= qpx_m_g(i_m_s+1)
                endif
                !u_o= (u2-u1)*(delx_m - delxj_m)/delxj_m + u2
                u_o= (u2-u1)*dx_s/delxj_m + u1
                if (ifc==1) then
                    q_g(i,j) = u_o
                elseif (ifc==2) then
                    qpx_g(i,j) = u_o
                endif
            enddo !* ifc=1,2
        enddo !* i=2, ncomp-1

    end subroutine mapback_qqpx

!*-----------------------------------------------------------*!






!**          Channel Geometry Calculation Tools             **!






!*-----------------------------------------------------------*!
    !+++---------------------------------------------------------------------------
    !+ Computation of area of various channel x-sections with depth as an argument
    !+ and without pre-specified chshp
    !+++---------------------------------------------------------------------------
    subroutine areacalc(hxs, Axs)

        implicit none
        double precision, intent(in) :: hxs
        double precision, intent(out) :: Axs
        double precision :: bwd, sslp, hbf, tw, twcc

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
        double precision, intent(in) ::  hxs, Axs
        double precision, intent(out) :: hydR
        double precision :: bxs, ssxs, twxs, twccxs
        double precision :: hbf

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
        double precision, intent(in) :: hxs, Axs, hydR
        double precision, intent(out) :: cnvey
        double precision :: bwd, sslp, Twd, TCCwd
        double precision :: subA, subP, hbf, TwCCi, K0, K1, K2

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
    !*--------------------------------------------------------------------------
    !*          Equivalent Manning's N for overbank flow in
    !*                                      trapz.main and rect. floodplain

    !*  Eq(4-35),p91,Chaudhry p95,RM1
    !*--------------------------------------------------------------------------
    double precision function emann_tmrf(hxs)
        implicit none
        !integer(KIND=i4b), intent(in) :: ic, jc
        double precision, intent(in) :: hxs
        double precision :: bwd, sslp, Twd, TCCwd
        double precision :: hbf, TwCCi, P0,P1,P2, tlP
        double precision ::  nom

        bwd= bo_g
        sslp= traps_g
        Twd= tw_g  !* top width of main channel
        TCCwd= twcc_g !* top width of compound channel
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
!*-----------------------------------------------------------*!






!**          Normal depth Calculation Tools             **!






!*-----------------------------------------------------------*!
!*-----------------------------------------------------------------------
!*              Bankfull depth-uniform bankfull discharge
!*                       for trapezoidal main channel
!*------------------------------------------------------------------------
    subroutine bfdep_ufQ(mxncomp_g, nrch_g, z_ar_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g,&
                        mann_ar_g, manncc_ar_g, dx_ar_g, so_ar_g, frnw_col, frnw_g)

        implicit none

        integer, intent(in) :: mxncomp_g, nrch_g, frnw_col
        double precision, dimension(mxncomp_g, nrch_g), intent(in) :: z_ar_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g
        double precision, dimension(mxncomp_g, nrch_g), intent(in) :: mann_ar_g, manncc_ar_g, dx_ar_g, so_ar_g
        integer, dimension(nrch_g, frnw_col), intent(in) :: frnw_g

        integer :: i, ncomp,j
        double precision :: hbf_tz, Axs, WP, hydR, ufQbf

        do j=1, nrch_g
            ncomp=frnw_g(j,1)
            do i=1, ncomp
                bo_g= bo_ar_g(i,j)
                traps_g= traps_ar_g(i,j)
                tw_g= tw_ar_g(i,j)
                mann_g= mann_ar_g(i,j)
                so_g= so_ar_g(i,j)

                !* bankfull depth
                hbf_tz= (tw_g - bo_g)/(2.0*traps_g)
                !* bankfull uniform discharge
                Axs=(bo_g  + traps_g*hbf_tz)*hbf_tz
                WP= bo_g + 2.0*hbf_tz*((1.0+(traps_g**2.0))**0.5)
                hydR=Axs/WP
                ufQbf= (1.0/mann_g)*Axs*(hydR**(2.0/3.0))*(so_g**0.5)

                bfuf_g(i,j,1)= hbf_tz
                bfuf_g(i,j,2)= ufQbf
            enddo
        enddo

    endsubroutine bfdep_ufQ
!*-----------------------------------------------------------------------
!*           Bisection method to compute normal depth using
!*            computed discharge for trapezoidal main channel, p4-1 & 5-1,RM5
!*------------------------------------------------------------------------
    subroutine bsec_nmdep(dsc, hbf, ufQbf, nmdep)

        implicit none

        double precision, intent(in) :: dsc, hbf, ufQbf
        double precision, intent(out) :: nmdep
        integer :: bs_flag, flrgm_flag
        double precision :: dscp, h, h0, h1, h2, h_dst0, h_dst, ufQ, ufQ1, ufQ2, dfQ0, dfQ1, dfQ2
        double precision :: cvrt, tolrt

!        open(unit=401, file="./output/bsec_nmdep.txt")

        tolrt=0.01
        dscp=abs(dsc)

        if (dscp.le.ufQbf) then
        !** Bisection for main channel
            h1=0.0
            h2=hbf
            flrgm_flag=1 !* indicate main channel flow
        else
        !** Bisection for flood plain
            !* first, find ufQ that is larger than dscp
            ufQ= ufQbf
            h= hbf
            do while (ufQ<dscp)
                h= h+hbf
                ufQ= ufQf(h, hbf)
            end do
            h1= hbf
            h2= h
            flrgm_flag=2 !* indicate flood plain flow
        endif

!            write(401,"(A10, 5A12, 2A20)") "inival:", "dsc", "hbf", "ufQbf",  "h1", "h2", "dfQ1", "dfQ2"
!                !* test only
!                if (flrgm_flag==1) then
!                !* main channel flow
!                    ufQ1= ufQm(h1)
!                    ufQ2= ufQm(h2)
!                elseif (flrgm_flag==2) then
!                !* floodplain flow
!                    ufQ1= ufQf(h1, hbf)
!                    ufQ2= ufQf(h2, hbf)
!                endif
!                dfQ1= ufQ1 -  dscp
!                dfQ2= ufQ2 -  dscp
!            write(401,"(A10, 5F12.6, 2F20.8)") "inival:", dsc, hbf, ufQbf,  h1, h2, dfQ1, dfQ2
!            write(401,"(A10, 5A12, 2A20, A12)") "bisec:", "dsc", "hbf", "ufQbf", "h1", "h2", "dfQ1", "dfQ2", "nmdep"

        h_dst0= h2-h1
        bs_flag=0
        do while (bs_flag==0)
            h0= 0.5*(h1+h2)
            if (flrgm_flag==1) then
                ufQ= ufQm(h0)
            elseif (flrgm_flag==2) then
                ufQ= ufQf(h0, hbf)
            endif

            dfQ0= ufQ -  dscp
            cvrt= abs(dfQ0)/ufQ

            if (cvrt.le.tolrt) then
                nmdep= h0
                bs_flag=1
                exit
            elseif (dfQ0.lt.0.0) then
                h1=h0
            elseif (dfQ0.gt.0.0) then
                h2=h0
            endif

            h_dst= h2-h1
            if (h_dst<0.05*h_dst0) then
            !* linear interpolation is used to wrap up the bisection
                if (flrgm_flag==1) then
                !* main channel flow
                    ufQ1= ufQm(h1)
                    ufQ2= ufQm(h2)
                elseif (flrgm_flag==2) then
                !* floodplain flow
                    ufQ1= ufQf(h1, hbf)
                    ufQ2= ufQf(h2, hbf)
                endif

                dfQ1= ufQ1 -  dscp
                dfQ2= ufQ2 -  dscp
                nmdep= -dfQ2*(h2-h1)/(dfQ2-dfQ1) + h2
                bs_flag=1
                !write(401,"(A10, 5F12.6, 2F20.8, F12.6)") "bisec:", dsc, hbf, ufQbf, h1, h2, dfQ1, dfQ2, nmdep
            end if
        enddo
    end subroutine bsec_nmdep
!*-----------------------------------------------------------------------
!*          Uniform main channel discharge for a given depth hxs
!*------------------------------------------------------------------------
    doubleprecision function ufQm(hxs)
        implicit none
        double precision, intent(in) :: hxs
        double precision :: Axs, WP, hydR

        !* trapezoidal channel inbank flow
        Axs=(bo_g + traps_g*hxs)*hxs
        WP= bo_g + 2.0*hxs*((1.0+(traps_g**2.0))**0.5)
        hydR=Axs/WP
        ufQm= (1.0/mann_g)*Axs*(hydR**(2.0/3.0))*(so_g**0.5)
    end function ufQm
!*-----------------------------------------------------------------------
!*          Uniform floodplain discharge for a given depth hxs
!*------------------------------------------------------------------------
    doubleprecision function ufQf(hxs, hbf)
        implicit none
        double precision, intent(in) :: hxs, hbf
        double precision :: Axs, WP, hydR, ufQ1, ufQ2, ufQ3

        !* subsection 1 and 3, p102,RM4
        Axs= (hxs-hbf)*(twcc_g-tw_g)/2.0
        WP= (twcc_g - tw_g)/2.0 + (hxs-hbf)
        hydR= Axs/WP
        ufQ1= (1.0/manncc_g)*Axs*(hydR**(2.0/3.0))*(so_g**0.5)
        ufQ3= ufQ1
        !* subsection 2, p102,RM4
        Axs= (bo_g + traps_g*hbf)*hbf + (hxs-hbf)*tw_g
        WP= bo_g + 2.0*hbf*((1.0 + traps_g**2.0)**0.5)
        hydR= Axs/WP
        ufQ2= (1.0/mann_g)*Axs*(hydR**(2.0/3.0))*(so_g**0.5)

        ufQf= ufQ1+ufQ2+ufQ3
    end function ufQf

endmodule mdiffv2