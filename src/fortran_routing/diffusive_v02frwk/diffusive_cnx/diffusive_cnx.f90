module diffusive
    !*----------------------------------------------------------------------------------------------------------
    !*      The CNX method for solving the diffusive wave equation is originally developed by Moussa R. and
    !*      et al. in apaper "Algorithms for solving the diffusive wave flood routing equation" and is
    !*      implemented on uniform space domain of each stream segment and integrated to a Python routing
    !*      framework by Dong Ha Kim, Adam Wlostowski, James Halgren, Jacob Hreha at the National Water Center. 
    !*      In the CNX, 
    !*      the PDE of diffusive equation is numerically solved by the Crank-Nicholson scheme with 2nd order 
    !*      solution in space. The derived numerical equations by C-N, which are a system of linear equations, 
    !*      are solved by resolution technique described in the paper. As CNX requires uniform space domain, 
    !*      each stream segment is sub-divided by a sub-segment length that is defined by the Courant condition.
    !*-----------------------------------------------------------------------------------------------------------
    implicit none

    double precision, parameter :: grav = 9.81
    double precision, parameter :: TOLERANCE = 1e-8
    integer :: nrch, mxncomp
    double precision :: dtini, t0, tfin, dt_ql, dt_ub, dt_db
    double precision, dimension(:,:), allocatable :: q_g, elv_g, qpx_g, q_pre_g, elv_pre_g
    double precision, dimension(:), allocatable :: tarr_ql, tarr_ub, tarr_db, varr_ql, varr_ub, varr_db
    integer, dimension(:,:), allocatable :: frnw_g
    double precision, dimension(:), allocatable :: qlatj_g
    double precision, dimension(:), allocatable :: qlatj_m_g
    double precision :: z_g, bo_g, traps_g, tw_g, twcc_g, so_g, mann_g, manncc_g
    double precision, dimension(:,:), allocatable :: adjz_ar_g, adjso_ar_g
    double precision, dimension(:,:), allocatable :: dx_ar, bo_ar, traps_ar, tw_ar, twcc_ar, mann_ar, manncc_ar
    integer, dimension(:,:), allocatable :: lcvQ_g, perQmap_g

    double precision, dimension(:,:), allocatable :: c_g, d_g
    integer :: perdim_g,  mxncomp_m_g
     double precision, dimension(:,:), allocatable :: bo_ar_m_g, traps_ar_m_g, tw_ar_m_g, twcc_ar_m_g
    double precision, dimension(:,:), allocatable :: mann_ar_m_g, manncc_ar_m_g, adjz_ar_m_g, dx_ar_m_g, adjso_ar_m_g
    double precision, dimension(:,:), allocatable :: q_m_g
    double precision, dimension(:,:), allocatable :: elv_m_g
    double precision, dimension(:,:), allocatable ::  qlatn_m_g
    double precision :: q_llm  !* minimum Q for numerical stability
    double precision, dimension(:,:), allocatable ::  mindepth_ns_ar_g
    double precision, dimension(:,:), allocatable :: mnq_drainql_g
    integer :: dim_nmLT_g
    double precision, dimension(:,:,:,:), allocatable :: normdepthLT_g, normdepthLT_m_g
    double precision, dimension(:,:), allocatable :: p, qq, r, w
    double precision, dimension(:,:), allocatable :: elv_diff_g
    double precision, dimension(:,:), allocatable :: dmy2d
    double precision, dimension(:), allocatable :: dmy1d
    integer, dimension(:,:), allocatable :: ncomp_mseg_g, oldncomp_mseg_g
    double precision, dimension(:,:), allocatable :: dx_mseg_g,   olddx_mseg_g
    double precision, dimension(:,:,:), allocatable :: q_mseg_g, elv_mseg_g, oldq_mseg_g, oldelv_mseg_g
    double precision, dimension(:,:,:), allocatable :: c_mseg_g, d_mseg_g
    double precision, dimension(:,:,:), allocatable :: p_mseg, qq_mseg, r_mseg, w_mseg
    integer :: nts_ub, mxncomp_mseg, ncomp_ghost
    double precision, dimension(:,:), allocatable :: qlat_tc, courant_comp
    double precision ::  C_llm, so_llm  !mn_chbtslp_c,
    double precision, dimension(:,:), allocatable :: q_mseg_btnode
    double precision, dimension(:), allocatable :: oldx, newx

contains
    !** Courant condition is applied to each segment so as to make uniform grid domain only within a segment not a reach to
    !** utilize the prismatic channel properties of a segment. Also, C and D are based on single-valued, stage-discharge
    !** relationship, P122,RM5.
    subroutine diffnw(timestep_ar_g, nts_ql_g, nts_ub_g, nts_db_g, ntss_ev_g, &
                      mxncomp_g, nrch_g, z_ar_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g, &
                      mann_ar_g, manncc_ar_g, so_ar_g, dx_ar_g, iniq, &
                      frnw_col, dfrnw_g, qlat_g, ubcd_g, dbcd_g, &
  		      paradim, para_ar_g, q_ev_g, elv_ev_g)
        implicit none

        integer, intent(in) :: mxncomp_g, nrch_g
        integer, intent(in) :: nts_ql_g, nts_ub_g, nts_db_g, ntss_ev_g
        integer, intent(in) :: frnw_col
        double precision, dimension(:), intent(in) :: timestep_ar_g(7)
        double precision, dimension(mxncomp_g, nrch_g), intent(in) :: z_ar_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g
        double precision, dimension(mxncomp_g, nrch_g), intent(in) :: mann_ar_g, manncc_ar_g, dx_ar_g, iniq
        double precision, dimension(nrch_g, frnw_col), intent(in) :: dfrnw_g
        double precision, dimension(nts_ql_g, mxncomp_g, nrch_g), intent(in) :: qlat_g
        double precision, dimension(nts_ub_g, nrch_g), intent(in) :: ubcd_g
        double precision, dimension(nts_db_g), intent(in) :: dbcd_g
	integer, intent(in) :: paradim
	double precision, dimension(paradim), intent(in) :: para_ar_g
        double precision, dimension(mxncomp_g, nrch_g), intent(inout) :: so_ar_g
        double precision, dimension(ntss_ev_g, mxncomp_g, nrch_g), intent(out) :: q_ev_g, elv_ev_g
        integer :: ncomp
        integer :: ts, n, i, j, k, ts_ev
        double precision ::  tc
        integer  ::  dsrchj
        double precision ::   saveInterval_min,  rmnd
        doubleprecision :: adjso, c_us, c_ds, so_us, so_ds, z_tw, hbf
        integer :: i2
        integer :: dimLT, ilt, ict, ict_1, ict_2, ict_3
        double precision :: step, qnorm, h_m, h_f
        integer :: twj, so_flag
        double precision :: Qxs, y_norm
        integer, dimension(:,:), allocatable :: pynw_g
        doubleprecision :: areai
        doubleprecision :: cfl, avec, sumc, vel_norm
        integer :: isubnode, subnode
        doubleprecision :: x1, x2
        doubleprecision :: newx_s, oldq1, oldq2, oldelv1, oldelv2  
	!*---------------------
	!* time step variables
	!*---------------------
        dtini= timestep_ar_g(1) 	!* simulation time step fixed throughout entire simulation time [sec]
        t0= timestep_ar_g(2) 		!* simulation start time [hr]
        tfin= timestep_ar_g(3) 		!* simulation end time [hr]
        !** not used: saveInterval= timestep_ar_g(4)  [sec]
        dt_ql= timestep_ar_g(5) 	!* lateral inflow data time step [sec]
        dt_ub= timestep_ar_g(6) 	!* upstream boundary discharge data time step [sec]
        dt_db= timestep_ar_g(7) 	!* downstream boundary stage data time step [sec]
	!*----------------------------
	!* sensitive model parameters
	!*----------------------------
	cfl = para_ar_g(1)  		!* upper limit of Courant number (default: 0.95) 
	C_llm = para_ar_g(2)  		!* lower limit of Celerity (default: 0.95) 
	q_llm = para_ar_g(3) 		!* lower limit of discharge (default: 0.002)
	so_llm = para_ar_g(4) 		!* lower limit of channel bed slope (default: 0.00001)
	mxncomp_mseg= int(para_ar_g(5)) !* the upper limit of number of uniformly distributed sub-nodes on each stream segment including ghost nodes (default:12)
	ncomp_ghost= int(para_ar_g(6))  !* the number of downstream sub-nodes added to existing sub-nodes on segment, used to provide proper discharge downstream boundary condition (default: 2)
        nrch=nrch_g
        mxncomp= mxncomp_g
        nts_ub = nts_ub_g

        allocate(dx_ar(mxncomp_g, nrch_g), bo_ar(mxncomp_g, nrch_g), traps_ar(mxncomp_g, nrch_g))
        allocate(tw_ar(mxncomp_g, nrch_g), twcc_ar(mxncomp_g, nrch_g))
        allocate(mann_ar(mxncomp_g, nrch_g), manncc_ar(mxncomp_g, nrch_g))
        dx_ar= dx_ar_g
        bo_ar= bo_ar_g
        traps_ar= traps_ar_g
        tw_ar= tw_ar_g
        twcc_ar= twcc_ar_g
        mann_ar= mann_ar_g
        manncc_ar= manncc_ar_g
        allocate(frnw_g(nrch_g,frnw_col))
        frnw_g=int(dfrnw_g)

        allocate(q_g(mxncomp_g, nrch_g), elv_g(mxncomp_g, nrch_g))
        allocate(q_pre_g(mxncomp_g, nrch_g), elv_pre_g(mxncomp_g, nrch_g))
        allocate(tarr_ql(nts_ql_g), varr_ql(nts_ql_g))
        allocate(tarr_ub(nts_ub_g), varr_ub(nts_ub_g))
        allocate(tarr_db(nts_db_g), varr_db(nts_db_g))
        allocate( ncomp_mseg_g(mxncomp_g, nrch_g), dx_mseg_g(mxncomp_g, nrch_g) )
        allocate( qlat_tc(mxncomp_g, nrch_g))
        allocate( courant_comp(mxncomp_g, nrch_g) )
        allocate( oldncomp_mseg_g(mxncomp_g, nrch_g), olddx_mseg_g(mxncomp_g, nrch_g) )
        allocate( q_mseg_g(mxncomp_mseg, mxncomp_g, nrch_g) )
        allocate(elv_mseg_g(mxncomp_mseg, mxncomp_g, nrch_g))
        allocate( oldq_mseg_g(mxncomp_mseg, mxncomp_g, nrch_g) )
        allocate( oldelv_mseg_g(mxncomp_mseg, mxncomp_g, nrch_g) )
        allocate( q_mseg_btnode(mxncomp_g, nrch_g) )
        allocate(elv_diff_g(mxncomp_g,nrch_g))

        !* time step series for lateral flow
        do n=1, nts_ql_g
            tarr_ql(n)= t0*60.0 + dt_ql*real(n-1,KIND(dt_ql))/60.0 !* [min]
        end do
        !* time step series for upstream boundary data
        do n=1, nts_ub_g
            tarr_ub(n)= t0*60.0 + dt_ub*real(n-1,KIND(dt_ub))/60.0 !* [min]
        enddo
        !* when measured data is used, time step series for downstream boundary data
!        do n=1, nts_db_g
!            tarr_db(n)= t0*60.0 + dt_db*real(n-1,KIND(dt_db))/60.0 !* [min]
!        enddo
        !**-----------------------------------------------------------------------------------*
        !*                  Adjust abnormally small channel bottom slope
        !*
        !* Definition of abnormal slope:
        !*  1. slopes less than a chosen lower limit (=so_llm)
        !* Repair method:
        !*  1. if less than lower limit, take average of slopes of adjacent segments
        !**-----------------------------------------------------------------------------------*
        allocate(adjso_ar_g(mxncomp_g,nrch_g))
        adjso_ar_g= so_ar_g
        !** NOTE:
        !* so_flag=1: activate the adjustment of bottom slope as below
        !* so_flag=0: Use slope as it is
        so_flag=0
        if (so_flag==1) then
            do j=1,nrch_g
                ncomp=frnw_g(j,1)
                do i=1, ncomp-1
                    if (so_ar_g(i,j).lt.so_llm) then
                        !* adjacent upstream segment's slope
                        so_us=0.0
                        c_us=0.0
                        i2=i-1
                        do while (i2.ge.1)
                            if (so_ar_g(i2,j).ge.so_llm) then
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
                            if (so_ar_g(i2,j).ge.so_llm) then
                                so_ds= so_ar_g(i2,j)
                                c_ds=1.0
                                exit
                            endif
                            i2=i2+1
                        end do
                        if (c_us+c_ds.gt.0.0) then
                            adjso= (so_us+so_ds)/(c_us+c_ds)
                            adjso_ar_g(i,j)= adjso
                        else
                            adjso_ar_g(i,j)= so_llm
                        endif
                    endif
                enddo
            enddo
        endif
        !++-----------------------------------------------------------------------------------+
        !+        adjust altitude one more time for counting for z_us = z_ds + so*dx
        !+
        !+        Refer to p48, RM5
        !++-----------------------------------------------------------------------------------+
        allocate(adjz_ar_g(mxncomp_g, nrch_g))
        do j=nrch_g, 1, -1
            if (frnw_g(j,2)<0.0) then
                twj= j !* TW reach index
                exit
            endif
        enddo
        ncomp=frnw_g(twj,1)
        z_tw= z_ar_g(ncomp, twj)
        do j=nrch_g, 1, -1
            ncomp=frnw_g(j,1)
            if (frnw_g(j,2)<0.0) then
            !* downstream boundary node at TW
                adjz_ar_g(ncomp, j)= z_tw
            else
            !* downstream boundary node at a junction
                dsrchj= frnw_g(j,2)    !* j index of the downstream reach
                adjz_ar_g(ncomp, j)= adjz_ar_g(1, dsrchj)
            endif
            !* For the following nodes within a reach before the bottom node:
            do i=ncomp-1,1,-1
                adjz_ar_g(i, j)= adjz_ar_g(i+1, j) + adjso_ar_g(i,j)*dx_ar_g(i,j)
            end do
        enddo
        !**----------------------------------------------------------------------*
        !*         Uniform flow lookup tables on original space domain
        !*
        !**----------------------------------------------------------------------*
        dim_nmLT_g=80
        dimLT=dim_nmLT_g
        allocate (normdepthLT_g(mxncomp_g,nrch_g,dimLT,2))
        do j=1, nrch_g
            do i=1, frnw_g(j,1)
                bo_g= bo_ar_g(i,j)
                traps_g= traps_ar_g(i,j)
                tw_g= tw_ar_g(i,j)
                twcc_g= twcc_ar_g(i,j)
                so_g= adjso_ar_g(i,j)
                mann_g= mann_ar_g(i,j)
                manncc_g= manncc_ar_g(i,j)

                hbf= (tw_g - bo_g)/(2.0*traps_g) !* bankfull depth
                !* main channel
                step= hbf/real(0.5*dimLT)
                h_m=0.0
                do ilt=1, int(0.5*dimLT)
                    call ufQ_tmrf(h_m, qnorm)
                    normdepthLT_g(i,j,ilt,1)= qnorm
                    normdepthLT_g(i,j,ilt,2)= h_m
                    h_m= h_m+ step
                enddo
                h_m= hbf
                call ufQ_tmrf(h_m, qnorm)
                ilt=int(0.5*dimLT)
                normdepthLT_g(i,j,ilt,1)= qnorm
                normdepthLT_g(i,j,ilt,2)= h_m
                !* flood plain
                step= 0.2*hbf
                h_f= hbf+step
                ict=0
                ict_1=int(0.25*(0.5*dimLT))
                ict_2=int(0.5*(0.5*dimLT))
                ict_3=int(0.75*(0.5*dimLT))
                do ilt=int(0.5*dimLT)+1, dimLT
                    call ufQ_tmrf(h_f, qnorm)
                    normdepthLT_g(i,j,ilt,1)= qnorm
                    normdepthLT_g(i,j,ilt,2)= h_f
                    ict=ict+1
                    if ((ict.gt.ict_1).and.(ict.le.ict_2)) then
                        step=0.5*hbf
                    elseif ((ict.gt.ict_2).and.(ict.le.ict_3)) then
                        step=1.0*hbf
                    elseif (ict.gt.ict_3) then
                        step=5.0*hbf
                    endif
                    h_f= h_f+ step
                enddo
            end do
        end do
        !**-----------------------------------------------------------------------------------*
        !*     INITIAL CONDITIONS of DISCHARGE on original space domain
        !*
        !**-----------------------------------------------------------------------------------*
        q_g= iniq
        !**-----------------------------------------------------------------------------------*
        !*                      Uniform domain of channel reach network
        !*
        !*  To apply CN method on uniform space domain, the number of new nodes are
        !*  computed by Courant condition. If computed number of nodes is less than three
         !**-----------------------------------------------------------------------------------*
        !* First, compute celerity using ck= (5/3)*V,V using Manning's Eq. <- Eq. (9.32),p346,Sturm.f
        allocate(c_g(mxncomp_g,nrch_g))
        do j=1, nrch_g
            do i=1, frnw_g(j,1)
                bo_g= bo_ar_g(i,j)
                traps_g= traps_ar_g(i,j)
                tw_g= tw_ar_g(i,j)
                twcc_g= twcc_ar_g(i,j)

                Qxs= q_g(i,j)
                y_norm= normdepth(i, j, Qxs)
                call areacalc(y_norm, areai)
                vel_norm= Qxs/areai
                c_g(i,j) = 5.0*vel_norm/3.0
                !* initial condition for water elev. on original space domain.
                elv_g(i,j) = y_norm + adjz_ar_g(i,j)
            enddo
            !* Second, assume celerities of nodes of the same reach ara all equal to each other.
            sumc=  sum(c_g(1:frnw_g(j,1),j))
            avec= sumc/real(frnw_g(j,1))
            do i=1, frnw_g(j,1)
                c_g(i,j)=max(avec, C_llm)
            enddo
        enddo

        call chgeo_seg_courant(cfl)

        deallocate(c_g)
        !**-----------------------------------------------------------------------------------*
        !*     INITIAL CONDITIONS of DISCHARGE on uniform space domain
        !*
        !**-----------------------------------------------------------------------------------*
        do j=1, nrch_g
            do i =1, frnw_g(j,1)-1
                do isubnode= 1, ncomp_mseg_g(i,j)
                    q_mseg_g(isubnode,i,j) = iniq(i,j)
                enddo
            enddo
        end do
        !**---------------------------------------------------------------------------------*
        !*                  INITIAL CONDITION of water depth

        !*      assuming depths of all reaches and nodes are equal to TW water depth.
        !**---------------------------------------------------------------------------------*
        do j=1, nrch_g
            do i =1, frnw_g(j,1)-1
                do isubnode= 1, ncomp_mseg_g(i,j)
                    Qxs= q_mseg_g(isubnode,i,j)
                    y_norm= normdepth(i, j, Qxs)
                    elv_mseg_g(isubnode,i,j) = y_norm +adjz_ar_g(i,j)
                end do
            enddo
        enddo
        !**---------------------------------------------------------------------------------*
        !*       minimum water depth responding to min.Q for numerical stability
        !*
        !**---------------------------------------------------------------------------------*
        allocate(mindepth_ns_ar_g(mxncomp_g, nrch_g))
        do j=1, nrch_g
            do i=1, frnw_g(j,1)
                mindepth_ns_ar_g(i,j)= normdepth(i, j, q_llm)
            end do
        end do

	ts=1        
	do j=1, nrch_g
            do i=1,frnw_g(j,1)
                q_ev_g(ts, i, j)=q_g(i,j)
                elv_ev_g(ts, i, j)= z_ar_g(i,j)+ elv_g(i,j)-adjz_ar_g(i,j)
            enddo
        end do

        allocate( c_mseg_g(mxncomp_mseg, mxncomp_g, nrch_g) )
        allocate( d_mseg_g(mxncomp_mseg, mxncomp_g, nrch_g) )
        allocate( p_mseg(mxncomp_mseg, mxncomp_g, nrch_g), qq_mseg(mxncomp_mseg, mxncomp_g, nrch_g) )
        allocate( r_mseg(mxncomp_mseg, mxncomp_g, nrch_g), w_mseg(mxncomp_mseg, mxncomp_g, nrch_g) )
        allocate( c_g(mxncomp_g,nrch_g), d_g(mxncomp_g,nrch_g) ) !* C and D variables of per matrix equation

        tc = t0*60.0  !* t0 is in hour. tc is in minutes

        do while (tc .lt. tfin*60.0)
            !**---------------------------------------------------------------------------------------------
            !*   discharge is computed by solving CNX with resolution technique on uniform grid, p91, RM5
            !**---------------------------------------------------------------------------------------------
            !* estimate lateral flow at current time n (=tc) for all segments of each reach.
            !* interpolated lateral flow for uniform space domain, p80,RM5
            qlat_tc = 0.0
            do j=1,nrch_g
                do i=1,frnw_g(j,1)-1
                    do n=1,nts_ql_g
                        varr_ql(n)= qlat_g(n,i,j) !* qlat_g(n,i,j) in unit of m2/sec
                    enddo
                    qlat_tc(i,j)= intp_y(nts_ql_g, tarr_ql, varr_ql, tc)
                enddo
            enddo

            call cnx_main(tc, ubcd_g) !*return q_m_g(i,j) for q at next time step.

            !* compute actual value of Courant number
            do j=1, nrch_g
                do i=1, frnw_g(j,1)-1
                    courant_comp(i,j) = c_g(i,j)*dtini/dx_mseg_g(i,j)
                enddo
            end do

            !**--------------------------------------------------------------------------------------------
            !* Map back q_m_g on uniform space domain to q_g on original domain  for elevation computation
            !**--------------------------------------------------------------------------------------------
            do j=1, nrch_g
                ncomp= frnw_g(j,1)
                i=1
                q_g(i,j)= q_mseg_g(1,i,j)
                !* option 1: take average of qs of sub-nodes between i-1 and i to compute q at i
!                do i=2,ncomp
!                    subnode_last= ncomp_mseg_g(i-1,j) - ncomp_ghost
!                    sumq= sum(q_mseg_g(1:subnode_last, i-1, j))
!                    aveq= sumq/real(subnode_last)
!                    q_g(i,j)= aveq
!                end do
!                !* make sure mass conserved at junction boundary
!                if (frnw_g(j,3)>0) then
!                    qjt= 0.0
!                    nusrch= frnw_g(j,3) !* then number of upstream reaches
!                    do rch=1, nusrch
!                        usrchj= frnw_g(j,3+rch) !* js corresponding to upstream reaches
!                        ncomp_usrchj= frnw_g(usrchj,1)
!                        qjt= qjt + q_g(ncomp_usrchj, usrchj)
!                    enddo
!                    q_g(1,j)= qjt
!                endif
                !* option 2: take instantaneous value of q at end sub-node of a segment b/t i-1 and i to compute q at i
                do i=2, ncomp-1
                    q_g(i,j)= q_mseg_g(1,i,j)
                enddo
                q_g(ncomp,j)= q_mseg_btnode(ncomp-1,j)
            enddo
            !*---------------------------------------------------------------------------------------
            !*                                       WATER DEPTH

            !* Compute y(=elv) or water depth at ts+1 from downstream reach to upstream reach while
            !* computing from downstream node to upstream node for each reach.
            !*---------------------------------------------------------------------------------------
            do j=nrch_g, 1, -1
                ncomp= frnw_g(j,1)
                !* downstream boundary condition for elv at ts+1
                if (frnw_g(j,2)<0.0) then
                !* downstream boundary node at TW
                !* 1. measured water depth data at TW node
!                    do n=1,nts_db_g
!                        varr_db(n)= dbcd_g(n) + adjz_ar_g(ncomp,j) !* when dbcd is water stage, channel bottom elev is added.
!                    enddo
!                    tf0= tc+dtini/60.0
!                    elv_g(ncomp,j)= intp_y(nts_db_g, tarr_db, varr_db, tf0)
                !* 2. use normal depth at  TW node
                    ncomp=frnw_g(j,1)
                    Qxs= q_g(ncomp,j)
                    y_norm= normdepth(ncomp, j, Qxs)
                    elv_g(ncomp,j) = y_norm + adjz_ar_g(ncomp,j)
                else
                !* downstream end node of a reach at a junction
                    !* 1. equal water depth at a juction
                    dsrchj= frnw_g(j,2)    !* j index for downstream reach
                    elv_g(ncomp,j)= elv_g(1,dsrchj)
                endif

                call elv_calc(ncomp, j) !* -> elv_g
            enddo
            !* compute celerity(ts+1) and diffusivity(ts+1) on the current sub-node configuration
            call C_D_mseg

            oldncomp_mseg_g = ncomp_mseg_g
            olddx_mseg_g = dx_mseg_g
            oldq_mseg_g = q_mseg_g

            call chgeo_seg_courant(cfl)

            !* compute q_mseg(t+1) and elv_mseg(t+1) on the new(ts+1) sub-node configuration by linear interpolation
            !* q_mseg needs to be computed upto ghost nodes, which is used in resolution method.
            do j=1, nrch_g
                do i=1, frnw_g(j,1)-1
                    allocate(oldx(oldncomp_mseg_g(i,j)))
                    do isubnode=1, oldncomp_mseg_g(i,j)
                        oldx(isubnode)= olddx_mseg_g(i,j)*real(isubnode-1)
                    enddo
                    !* find subnodes of oldx where new subnode location is inbetween
                    do isubnode=1, ncomp_mseg_g(i,j)
                        newx_s = dx_mseg_g(i,j)*real(isubnode-1)
                        subnode =locate(oldx,newx_s) !* oldx(subnode)<newx_s<oldx(subnode+1)

                        if (subnode==oldncomp_mseg_g(i,j)) then
                            x1= oldx(subnode-1); oldq1= oldq_mseg_g(subnode-1,i,j)
                            x2= oldx(subnode); oldq2= oldq_mseg_g(subnode,i,j)
                        else
                            x1= oldx(subnode); oldq1= oldq_mseg_g(subnode,i,j)
                            x2= oldx(subnode+1); oldq2= oldq_mseg_g(subnode+1,i,j) 
                        endif
                        q_mseg_g(isubnode,i,j) = LInterpol(x1, oldq1, x2, oldq2, newx_s)
                    enddo
                    deallocate(oldx)
                enddo
            enddo

            tc= tc + dtini/60.0 !* [min]
            ts= ts+1

            do j=1, nrch_g
                ncomp= frnw_g(j,1)
                do i=1, ncomp
                    q_ev_g(ts, i, j)=q_g(i,j)
                    elv_ev_g(ts, i, j)= z_ar_g(i,j)+ elv_g(i,j)-adjz_ar_g(i,j)
                enddo
            enddo
        enddo !* end of simulation time loop

        deallocate(dx_ar, bo_ar, traps_ar, tw_ar, twcc_ar, mann_ar, manncc_ar)
        deallocate(frnw_g)
        deallocate(q_g, elv_g, q_pre_g, elv_pre_g)
        deallocate(tarr_ql, varr_ql, tarr_ub, varr_ub, tarr_db, varr_db)
        deallocate(adjso_ar_g, adjz_ar_g)
        deallocate(normdepthLT_g)
        deallocate(c_g, d_g)
        deallocate(mindepth_ns_ar_g)
        deallocate( p_mseg, qq_mseg, r_mseg, w_mseg)

    endsubroutine diffnw
    !**------------------------------------------------------------------------------
    !*         CNX and resolution technique on uniform grid, p106~109,RM5

    !**------------------------------------------------------------------------------
    subroutine cnx_main(tc, ubcd_g)
        implicit none
        doubleprecision, intent(in) :: tc
        double precision, dimension(nts_ub, nrch), intent(in) :: ubcd_g
        !* initialize
        p_mseg = 0.0
        qq_mseg = 0.0
        r_mseg = 0.0
        w_mseg = 0.0
        !* Compute C and D at (i,j) and n, p27_10~11
        call C_D_mseg
        !* Compute p, q, r, and w coefficients for entire network at time n
        call p_q_r_w_mseg
        !* Compute Q using resolution technique
        call Q_mseg_resoltech(tc, ubcd_g)

     endsubroutine cnx_main
    !**------------------------------------------------------------------------------
    !*         Celerity and Diffusivity assuming uniform flow along a segment

    !*                              p122, RM5
    !**------------------------------------------------------------------------------
    subroutine C_D_mseg
        implicit none
        integer :: i, j
        double precision ::  avec, aved, sumc, sumd
        double precision ::  hxs, Axs, Bxs, Qxs, chbed_slp, emann 
        integer :: isubnode, subnode_last
        double precision :: c_nom, c_denom
 
        do j=1, nrch
            do i=1, frnw_g(j,1)-1
                bo_g= bo_ar(i,j)
                traps_g= traps_ar(i,j)
                tw_g= tw_ar(i,j)
                twcc_g= twcc_ar(i,j)
                mann_g= mann_ar(i,j)
                manncc_g= manncc_ar(i,j)
                chbed_slp = max(so_llm, adjso_ar_g(i,j))

                subnode_last = ncomp_mseg_g(i,j) -1
                do isubnode= 1, subnode_last
                    !* compute area and top width at each sub-node of a segment                   
                    Qxs= q_mseg_g(isubnode,i,j)
		    hxs= normdepth(i, j, Qxs)
		    hxs= max(hxs, mindepth_ns_ar_g(i,j))
                    call areacalc(hxs, Axs)
                    Bxs= topwidth(hxs)                    
                    !* celerity at sub-node isubnode at i segment of j reach
                    !* option 1
                    c_mseg_g(isubnode,i,j) = (5.0/3.0)*Qxs/Axs
                    !* option 2: Celerity by a paper of Litrico et al.
!                    emann= emann_tmrf(hxs) !* equivalent manning's N
!                    c_nom = (chbed_slp**0.3)*(q_mseg_g(isubnode,i,j)**0.4)
!                    c_denom = (Bxs**0.4)*(emann**0.6)
!                    c_mseg_g(isubnode,i,j) = (5.0/3.0)*c_nom/c_denom
                    !* diffusivity at sub-node isubnode at i segment of j reach
                    d_mseg_g(isubnode,i,j) = Qxs/(2.0*Bxs*chbed_slp)
                enddo
                !* option 1
                !* use diffusivity and celerity values averaged along all the sub-nodes of a segment
                sumc=  sum(c_mseg_g(1:subnode_last, i, j))
                avec= sumc/real(subnode_last)
                c_g(i,j) = max(avec, C_llm)
                sumd= sum(d_mseg_g(1:subnode_last, i, j))
                aved= sumd/real(subnode_last)
                d_g(i,j) = aved
            enddo
        enddo
    endsubroutine C_D_mseg
    !**-------------------------------------------------------------------*
    !*       p, q, r, w coefficients of CNX along sub-nodes of a segment

    !*                          p90, RM5
    !**-------------------------------------------------------------------*
     subroutine p_q_r_w_mseg
        implicit none
        integer :: i, j
        double precision :: dxst, hhi, ggi, ppr, qpr, rpr, avec, aved
        integer :: isubnode, subnode_last

        do j=1, nrch
            do i=1, frnw_g(j,1)-1
                dxst = dx_mseg_g(i,j)
                subnode_last = ncomp_mseg_g(i,j) - ncomp_ghost
                do isubnode= 2, ncomp_mseg_g(i,j)-1
                    avec = c_g(i,j)
                    aved = d_g(i,j)

                    hhi= avec*dtini/dxst
                    ggi= aved*dtini/(dxst**2.0)

                    p_mseg(isubnode,i,j)= -hhi/4.0 - ggi/2.0
                    qq_mseg(isubnode,i,j)= 1.0 + ggi
                    r_mseg(isubnode,i,j)= hhi/4.0 - ggi/2.0

                    ppr= hhi/4.0 + ggi/2.0
                    qpr= 1.0 - ggi
                    rpr= -hhi/4.0 + ggi/2.0

                    w_mseg(isubnode, i,j)= ppr*q_mseg_g(isubnode-1,i,j) + qpr*q_mseg_g(isubnode,i,j) &
                            + rpr*q_mseg_g(isubnode+1,i,j) + avec*qlat_tc(i,j)*dtini
                enddo
            enddo
        enddo
    endsubroutine p_q_r_w_mseg
    !**-----------------------------------------------------------------------------------------------------*
    !*              Q computation using resolution technique of CNX, p106,RM5 (Q_ncomp(n+1) = Q_ncomp-1(n+1)

    !* source: 1. Moussa at el. "Algorithms for solving the diffusive wave flood routing equation",
    !*              Hydrological Processes. vol 10, 105-123 (1996)
    !**------------------------------------------------------------------------------------------------------*
     subroutine Q_mseg_resoltech(tc, ubcd_g)
        implicit none
        doubleprecision, intent(in) :: tc
        double precision, dimension(nts_ub, nrch), intent(in) :: ubcd_g
        integer :: i, j, nusrch, rch, usrchj
        double precision :: dxst
        double precision :: qjt
        double precision, dimension(:), allocatable :: x, y
        integer :: ncomp_ms
        integer :: isubnode, n, subnode_last, subnode
        double precision :: tf0,totalx, x1, x2, q1_mseg, q2_mseg

        do j=1, nrch
            do i =1, frnw_g(j,1)-1
                dxst = dx_mseg_g(i,j)
                ncomp_ms = ncomp_mseg_g(i,j)
                if (i==1) then
                    if (frnw_g(j,3)==0) then !* frnw_g(j,3) shows the number of upstream reaches.
                    !* head water reach
                        do n=1,nts_ub
                            varr_ub(n)= ubcd_g(n,j)
                        enddo
                        tf0= tc+dtini/60.0
                        q_mseg_g(1,i,j)= intp_y(nts_ub, tarr_ub, varr_ub, tf0) !* tarr_ub in min.
                    else
                    !* at a junction
                        qjt= 0.0
                        nusrch= frnw_g(j,3) !* then number of upstream reaches
                        do rch=1, nusrch
                            usrchj= frnw_g(j,3+rch) !* js corresponding to upstream reaches
                            qjt = qjt + q_mseg_btnode(frnw_g(usrchj,1)-1, usrchj)
                        enddo
                         q_mseg_g(1,i,j)= qjt
                    endif
                    !* check against min.q for stability
                    q_mseg_g(1,i,j)= max(q_mseg_g(1,i,j), q_llm)
                endif
                !* Q at the top sub-node of from 2nd segment to onward is equal to Q at the last sub-node of its own upstream segment.
                if (i.ge.2) then
                    q_mseg_g(1,i,j) = q_mseg_btnode(i-1,j)
                endif
               !* x and y computation in backward-sweep
                allocate(x(ncomp_ms), y(ncomp_ms))
                x(ncomp_ms)= 1.0
                y(ncomp_ms)= 0.0
                do isubnode=ncomp_ms-1, 2, -1
                    x(isubnode)= -p_mseg(isubnode,i,j)/(qq_mseg(isubnode, i,j) + r_mseg(isubnode, i,j)*x(isubnode+1))
                    y(isubnode)= (w_mseg(isubnode, i,j)- r_mseg(isubnode,i,j)*y(isubnode+1))/(qq_mseg(isubnode,i,j)+ &
                                    r_mseg(isubnode,i,j)*x(isubnode+1))
                end do
                !* Q computation at time n+1 in forward-sweep
                do isubnode=2, ncomp_ms
                    q_mseg_g(isubnode,i,j)= x(isubnode)*q_mseg_g(isubnode-1,i,j) + y(isubnode)
                    q_mseg_g(isubnode,i,j)= max(q_mseg_g(isubnode,i,j), q_llm)
                enddo
                deallocate(x, y)

                totalx= dx_ar(i,j)
                allocate (oldx(ncomp_ms))
                do isubnode= 1, ncomp_ms
                    oldx(isubnode) = dxst*real(isubnode-1)
                end do
                subnode =locate(oldx,totalx) !* oldx(subnode)<totalx<oldx(subnode+1)
                if (subnode==ncomp_ms) then
                    x1= oldx(subnode-1); q1_mseg= q_mseg_g(subnode-1,i,j)
                    x2= oldx(subnode); q2_mseg= q_mseg_g(subnode,i,j)
                else
                    x1= oldx(subnode); q1_mseg= q_mseg_g(subnode,i,j)
                    x2= oldx(subnode+1); q2_mseg= q_mseg_g(subnode+1,i,j)
                endif
                q_mseg_btnode(i, j) = LInterpol(x1, q1_mseg, x2, q2_mseg, totalx)
                deallocate(oldx)
            enddo
        enddo
    endsubroutine Q_mseg_resoltech
    !*--------------------------------------
    !*              conveyance
    !*-------------------------------------
    doubleprecision function convey(i, j, depthi)
        implicit none
        integer, intent(in) :: i, j
        double precision, intent(in) :: depthi
        double precision :: areai, hydR, conv
        !* conveyance at (i,j,n)
        bo_g= bo_ar(i,j)
        traps_g= traps_ar(i,j)
        tw_g= tw_ar(i,j)
        twcc_g= twcc_ar(i,j)
        mann_g= mann_ar(i,j)
        manncc_g= manncc_ar(i,j)

        call areacalc(depthi, areai)
        call hydRcalc(depthi, areai, hydR)
        call Kcalc(depthi, areai, hydR, conv)
        convey= conv
    endfunction convey
    !**--------------------------------------------------------------
    !*          Computation of top width of water surface
    !**--------------------------------------------------------------
    doubleprecision function topwidth(hxs)
        implicit none
        doubleprecision, intent(in) :: hxs
        doubleprecision :: bwd, sslp, hbf, tw, twcc

        bwd= bo_g !*bottom width
        sslp= traps_g !*trapezoidal
        tw= tw_g
        twcc=twcc_g
        if (sslp==0.0) then
        !* rectangular channel
            topwidth= bwd
        else
            hbf= (tw - bwd)/(2.0*sslp) !* bankfull depth
            if (hxs.le.hbf) then
                !* trapezoidal channel inbank flow
                topwidth= bwd + 2.0*sslp*hxs
            else
                !*overbank flow on rect. floodplains
                topwidth= twcc
            end if
        endif
    endfunction topwidth
    !*--------------------------------------------------------------------------
    !*          Equivalent Manning's N for overbank flow in
    !*                                      trapz.main and rect. floodplain

    !*  Eq(4-35),p91,Chaudhry p95,RM1
    !*--------------------------------------------------------------------------
    doubleprecision function emann_tmrf(hxs)
        implicit none
        doubleprecision, intent(in) :: hxs
        doubleprecision :: bwd, sslp, Twd, TCCwd, mann, manncc
        doubleprecision :: hbf, TwCCi, P0,P1,P2, tlP
        doubleprecision ::  nom

        bwd= bo_g
        sslp= traps_g
        Twd= tw_g
        TCCwd= twcc_g
        mann= mann_g
        manncc= manncc_g
        !* bankfull depth
        hbf=(Twd - bwd)/(2.0*sslp)
        if (hxs.gt.hbf) then
            !*P0,P1,P2 for overbank flow p84-2-1,RM1
            P0 = bwd + 2.0*hbf*(1.0 + sslp**2.0)**0.5
            TwCCi= (TCCwd-Twd)/2.0
            P1=TwCCi +  (hxs - hbf)
            P2=TwCCi +  (hxs - hbf)

            tlP= P0+P1+P2
            !* equivalent Manning's N for compound channel, Eq(4-35),p91,Chaudhry
            nom= P0*mann**(3.0/2.0) + P1*manncc**(3.0/2.0) + P2*manncc**(3.0/2.0)
            emann_tmrf= (nom/tlP)**(2.0/3.0)
        else
            !* still in main channel
            emann_tmrf= mann
        endif
    end function emann_tmrf
    !*----------------------------------------------------------------------------------------
    !*      Compute water elevation from downstream to upstream based on normal depth

    !*-----------------------------------------------------------------------------------------
    subroutine elv_calc(ncomp, j)
        implicit none
        integer,intent(in) :: ncomp, j
        integer :: i
        doubleprecision :: qi, yi_tp1, sf, convk, dip1

        !* compute water elevation at the next time step
        do i=ncomp-1,1,-1
            qi= q_g(i,j)
            yi_tp1= normdepth(i, j, qi)
            !* make sure water depth larger than arbitrarily chosen minimum positive water depth.
            yi_tp1= max(yi_tp1, mindepth_ns_ar_g(i,j))
            elv_g(i,j)= yi_tp1 + adjz_ar_g(i,j)
        end do
    end subroutine elv_calc
    !**----------------------------------------------------------------
    !*      compute normal depth using Bisection method

    !*-----------------------------------------------------------------
    doubleprecision function normdepth(i, j, Qxs)
        implicit none
        integer, intent(in) :: i, j
        double precision, intent(in) :: Qxs
        double precision, dimension(:), allocatable :: xarr, yarr
        integer :: ilt, dimLT

        dimLT= dim_nmLT_g
        allocate(xarr(dimLT), yarr(dimLT))
        do ilt=1, dimLT
            xarr(ilt)= normdepthLT_g(i,j,ilt,1) !* uniform discharge
            yarr(ilt)= normdepthLT_g(i,j,ilt,2) !* normal depth
        enddo
        normdepth= intp_y(dimLT, xarr, yarr, Qxs)
        deallocate(xarr,yarr)
    endfunction normdepth
    !*--------------------------------------------
    !           Interpolate y for given x
    !
    !*--------------------------------------------
    doubleprecision function intp_y(nrow, xarr, yarr, x)
        implicit none
        integer, intent(in) :: nrow
        doubleprecision, dimension(nrow), intent(in) :: xarr, yarr
        doubleprecision, intent(in) :: x
        integer :: irow
        doubleprecision :: x1, y1, x2, y2, y

        irow= locate(xarr, x)
        if (irow.eq.0) irow= 1
        if (irow.eq.nrow) irow= nrow-1
        x1= xarr(irow); y1= yarr(irow)
        x2= xarr(irow+1); y2= yarr(irow+1)
        y= LInterpol(x1,y1,x2,y2,x)
        intp_y = y
    end function intp_y
    !*------------------------------------------------------------------------------------
    !                       Locate function in f90, p.1045,NR f90
    !
    !   klo=max(min(locate(xa,x),n-1),1) In the Fortran 77 version of splint,
    !   there is in-line code to find the location in the table by bisection. Here
    !   we prefer an explicit call to locate, which performs the bisection. On
    !   some massively multiprocessor (MMP) machines, one might substitute a different,
    !   more parallel algorithm (see next note).
    !*------------------------------------------------------------------------------------
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
    doubleprecision function LInterpol(x1,y1,x2,y2,x)
        implicit none
        doubleprecision, intent(in) :: x1, y1, x2, y2, x
        !* interpolate y for the given x
        LInterpol= (y2-y1)/(x2-x1)*(x-x1)+y1
    end function LInterpol
    !*----------------------------------------------------------------------
    !*     Uniform flow discharge for trapz. main and rect. floodplain
    !*     p.102, RM4
    !*----------------------------------------------------------------------
    subroutine ufQ_tmrf(hxs, ufQ)
        implicit none
        double precision, intent(in) :: hxs
        double precision, intent(out) :: ufQ
        double precision :: Axs, hydR, WP, Twd, TCCwd, So
        double precision :: bwd, sslp, hbf, ufQ1, ufQ2, ufQ3
        bwd= bo_g
        sslp= traps_g
        Twd= tw_g       !* top width of main channel
        TCCwd= twcc_g   !* top width of compound channel
        So= so_g
        hbf= (Twd - bwd)/(2.0*sslp) !* bankfull depth
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
    !**------------------------------------------------------------------------------
    !*  Computation of area of various channel x-sections with depth as an argument
    !*
    !**------------------------------------------------------------------------------
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
            !* rectangular channel
            Axs=bwd*hxs
        else
            hbf= (tw - bwd)/(2.0*sslp) !* bankfull depth
            if (hxs.le.hbf) then
                !* trapezoidal channel inbank flow
                Axs=(bwd + sslp*hxs)*hxs
            else
                !*overbank flow on rect. floodplains
                Axs=(bwd + sslp*hbf)*hbf + twcc*(hxs-hbf)
            end if
        endif
    end subroutine areacalc
    !**-----------------------------------------------------
    !*      Computation of hydraulic radius R (=A/P)
    !**-----------------------------------------------------
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
    !**-----------------------------------------------------
    !* Computation of conveyance K (=1/N * A * R^(2/3))
    !**-----------------------------------------------------
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
    endsubroutine Kcalc
    !*----------------------------------------------------------------------------------------------
    !* construct uniform hydrofabric for each segment not each reach according to Courant condition
    !*----------------------------------------------------------------------------------------------
    subroutine chgeo_seg_courant(courantN)
        implicit none
        double precision, intent(in) :: courantN
        integer :: i, j, ncomp, isubnode
        doubleprecision :: delx_cr, sumdx,  totalx_mseg, delx_cr_final, incrdelx_cr

        do j=1,nrch
            !* sub-segment length defined by Courant condition
            ncomp= frnw_g(j,1)
            do i=1, ncomp-1
                delx_cr = c_g(i,j)*dtini/courantN  !* assumed celerity is same for all nodes of the same reach.
                sumdx= dx_ar(i,j)
                !* when delx_cr is too small, adjust it according to mxncomp_mseg
                incrdelx_cr= 0.5*delx_cr
                delx_cr_final= delx_cr
                totalx_mseg = delx_cr_final*(mxncomp_mseg-ncomp_ghost-1)
                do while (totalx_mseg.lt.sumdx)
                    delx_cr_final = delx_cr_final + incrdelx_cr
                    totalx_mseg = delx_cr_final*(mxncomp_mseg-ncomp_ghost-1)
                end do
                delx_cr = delx_cr_final
                !* new sub-node configuration is constructed. p129,RM5
                totalx_mseg=0.0
                isubnode=0
                do while (totalx_mseg.le.sumdx)
                    isubnode= isubnode + 1
                    totalx_mseg= delx_cr*real(isubnode-1)
                enddo
                ncomp_mseg_g(i,j) = isubnode + ncomp_ghost
                dx_mseg_g(i,j)= delx_cr
            enddo
        enddo
    endsubroutine chgeo_seg_courant
end module diffusive
