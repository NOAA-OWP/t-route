!            FINITE DIFFERENCE METHOD
!
!  A program for one dimensional flow in open channel
!
program mesh

    use var
    use subtools
    use dimlessN
    use tstep
    use mesh_kernel
    use lookupTB

    implicit none

    ! Local storage
    integer :: i, j, k, ppn, qqn, n, ntim, igate, pp, boundaryFileMaxEntry, saveFrequency
    integer :: linknb_ds, linknb_us
    real :: qnp1_ds, qnp1_us, qsum, y_ds
    real :: cour, da, dq, x, saveInterval, width
    real :: qn, xt, maxCourant, dtini_given, nodenb, linknb
    real :: frds, areasum, yk_ncomp, yav, areak_ncomp, areav, sumOldQ, currentQ, area_ds
    real :: arean, areac, hyrdn, hyrdc, perimn, perimc, qcrit, s0ds, timesDepth, latFlowValue, latFlowValue2
    real :: t, r_interpol_time, tfin, t1, t2, t0 !t0 start time
    integer :: tableLength, timestep, kkk
    real :: area_0, width_0, errorY, hydR_0, q_sk_multi, sumCelerity
    real :: r_interpo_nn
    real :: maxCelDx
    real, dimension(:), allocatable :: dmyv, dmyv1, dmyv2
    real, dimension(:), allocatable ::tarr_ql, varr_ql, tarr_ub, varr_ub, tarr_db, varr_db
    integer :: frj, iseg, i1
    integer :: num_points, totalChannels
    real, allocatable :: leftBank(:,:), rightBank(:,:)
    real, allocatable :: skLeft(:,:), skMain(:,:), skRight(:,:)
    real :: dmy1, dmy2
    integer :: ndata, idmy1
    integer(8), allocatable :: pynw(:,:)

    call cpu_time( t1 )
!**
!* START: v02 routing framework input data
!**
    nrch_g=735  !29   !735             !* # of total reaches TW:933020089  USGS:02086624
    dtini_g= 10.0         !* initial time interval in [sec]
    t0_g= 0.0             !* simulation starting time in [hr]
    tfin_g= 718 !408 !718.0       !*simulation ending time in [hr] TW:933020161 USGS:02086849
    tfin_qlat_g=718.0   !*added as qlat data time period is generally not the same as dsbd data. TW:933020161 USGS:02086849.
    saveInterval_g=dtini_g    !* output publishing time interval in [sec]
    cfl_g= 0.99            !* max. allowable Courant number
    dt_ql_g= 3600.0  !* instead of dtinput #time interval of qlateral input data from wrf-hydro [sec]
    dt_ub_g= dt_ql_g  !* time interval of input data for upstream boundary condition  [sec]
    dt_db_g=  900.0 !* time interval of input data for downstream boundary condition  [sec]
    saveInterval_ev_g= 1800 !7200 !saveInterval_g !for evaluation later [sec]
    nl_ubcd_ar_g=266749 !10785  !266749  !* # of total rows of ubcd_ar.txt TW:933020089  USGS:02086624

    ntss_ev_g= int((tfin_g - t0_g)*3600.0/saveInterval_ev_g, KIND(ntss_ev_g))+4 !*# of timesteps for publishing final outputs
    nts_ql_g= int( (tfin_qlat_g - t0_g)*3600.0/dt_ql_g+1, KIND(nts_ql_g) ) !* # of timesteps of lateral inflow input data
    nts_ub_g= nts_ql_g  !* # of timesteps of upstream boundary input data
    nts_db_g= int( (tfin_g - t0_g)*3600.0/dt_db_g+1, KIND(nts_db_g) ) !* # of timesteps of downstream boundary input data

    timesdepth_g=4.0 !* water depth multiplier used in readXsection
    nel_g=501         !* # of sub-Xsection used to construct hydraulic lookup table.    ! change Nazmul 20210601

    open(unit=21, file="./input/frnw_ar.txt")
    open(unit=22, file="./input/z_bo_traps_tw_twcc_m_mcc_so_dx.txt")
    open(unit=22, file="./input/z_ar.txt")
    open(unit=23, file="./input/bo_ar.txt")
    open(unit=24, file="./input/traps_ar.txt")
    open(unit=25, file="./input/tw_ar.txt")
    open(unit=26, file="./input/twcc_ar.txt")
    open(unit=27, file="./input/mann_ar.txt")
    open(unit=28, file="./input/manncc_ar.txt")
    open(unit=29, file="./input/so_ar.txt")
    open(unit=30, file="./input/dx_ar.txt")
    open(unit=31, file="./input/qlat_ar.txt")
    open(unit=32, file="./input/ubcd_ar.txt")
    !open(unit=33, file="./input/dbcd_ar.txt")
    open(unit=33, file="./input/dbcd_ar_datumCorrected.txt")    ! change Nazmul 20210601    !! The previous d/s boundary was not datum corrected
    open(unit=34, file="./input/pynw_ar.txt")
    open(unit=101, file="./output/simulated discharge depth elev.txt")

    nlinks=nrch_g
    !* frnw contains information about channel connectivity.
    !* First of all, let's call j reach index, where j=1,..,# of total reaches, that corresponds to the ID of head link of a reach.
    !* The value of j is determined by the number of junctions downstream of a individual reach.  So, for example, when a reach has
    !* the maximum number of junctions in the downstream, then its j is 1. If there are multiple reaches with the max. value, then
    !* they are randomly assigned values between 1 and the number of the reaches.  When it is only one that has no junction downstream,
    !* its j value becomes # of total reaches.
    !* frnw(j,1) = # of nodes of reach j.  Here, j actually corresponds to the head node of the reach.
    !* frnw(j,2) = reach index corresponding to the adjacent downstream reach (or head node of the reach) of reach j.
    !* frnw(j,3)= # of reaches in the adjacent upstream of reach j
    !* frnw(j, 3+i) = reach indices corresponding to individual upstream reaches.
    !* For example, reach 1 and 2 together joins at reach 3. Reach 6 is right below reach 3. Reach 3 has eight nodes.
    !* The link ID of the top node of reach 3 is 880423.
    !* Then,
    !* frnw(3,1)=8
    !* frnw(3,2)=6
    !* frnw(3,3)=2
    !* frnw(3,4)=1
    !* frnw(3,5)=2
    allocate(frnw_g(nlinks,8), dfrnw_g(nlinks,8))
    do j=1, nlinks
        read(21,*) (dfrnw_g(j,i),i=1,8)
    enddo
    frnw_g=dfrnw_g
    allocate(dmyv(nlinks))
    dmyv=frnw_g(:,1)
    mxncomp_g= maxval(dmyv)
    !* pynw maps j values back into link ID of head node of a reach
    !* So, from the previous example, for reach 3 pynw(3,1)=3 and pynw(3,2)=880423
    allocate(pynw(nlinks,2))
    do j=1, nlinks
        read(34,*) pynw(j,1), pynw(j,2)
    enddo

    allocate( z_ar_g(mxncomp_g, nlinks), bo_ar_g(mxncomp_g, nlinks), traps_ar_g(mxncomp_g, nlinks) )
    allocate( tw_ar_g(mxncomp_g, nlinks), twcc_ar_g(mxncomp_g, nlinks) )
    allocate( mann_ar_g(mxncomp_g, nlinks), manncc_ar_g(mxncomp_g, nlinks) )
    allocate( so_ar_g(mxncomp_g, nlinks), dx_ar_g(mxncomp_g, nlinks) )
    allocate( qlat_g(nts_ql_g, mxncomp_g, nlinks) )
    allocate(tarr_ql(nts_ql_g), varr_ql(nts_ql_g))
    allocate(tarr_ub(nts_ub_g), varr_ub(nts_ub_g))

    do i=1, mxncomp_g
        read(22,*) (z_ar_g(i,j),j=1,nlinks)
    enddo
    do i=1, mxncomp_g
        read(23,*) (bo_ar_g(i,j),j=1,nlinks)
    enddo
    do i=1, mxncomp_g
        read(24,*) (traps_ar_g(i,j),j=1,nlinks)
    enddo
    do i=1, mxncomp_g
        read(25,*) (tw_ar_g(i,j),j=1,nlinks)
    enddo
    do i=1, mxncomp_g
        read(26,*) (twcc_ar_g(i,j),j=1,nlinks)
    enddo
    do i=1, mxncomp_g
        read(27,*) (mann_ar_g(i,j),j=1,nlinks)
    enddo
    do i=1, mxncomp_g
        read(28,*) (manncc_ar_g(i,j),j=1,nlinks)
    enddo
    do i=1, mxncomp_g
        read(29,*) (so_ar_g(i,j),j=1,nlinks)
    enddo
    do i=1, mxncomp_g
        read(30,*) (dx_ar_g(i,j),j=1,nlinks)
    enddo
    do j=1, nrch_g
        ncomp=frnw_g(j,1)
        do i=1, ncomp
            do n=1, nts_ql_g
                read(31,*) frj, iseg, tarr_ql(n), qlat_g(n,i,j) !* qlat_g(n,i,j) in [m^2/sec]
            end do
        end do
    enddo
!**
!* END: v02 routing framework input data
!**
    dtini= dtini_g
    dtini_given= dtini
    lastKnownDiffuDT = dtini
    t0=t0_g
    tfin= tfin_g
	ntim = floor( (tfin - t0) / dtini * 3600)
    timesDepth= timesdepth_g
    nel= nel_g
    saveInterval= saveInterval_ev_g
    saveFrequency = saveInterval / dtini_given
    num_points= mxncomp_g
    totalChannels= nlinks
    allocate(area(num_points))
    ! change for unsteady flow
    allocate(bo(num_points,totalChannels))
    allocate(pere(num_points,totalChannels))
    allocate(dpda(num_points))
    allocate(areap(num_points,totalChannels))
    allocate(qp(num_points,totalChannels))
    allocate(z(num_points,totalChannels))
    allocate(dqp(num_points,totalChannels))
    allocate(dqc(num_points,totalChannels))
    allocate(dap(num_points,totalChannels))
    allocate(dac(num_points,totalChannels))
    allocate(depth(num_points))
    allocate(sk(num_points,totalChannels))
    allocate(co(num_points))
    allocate(dbdx(num_points))
    allocate(dt(num_points))
    allocate(dx(num_points,totalChannels))
    allocate(volRemain(num_points-1,totalChannels))
    allocate(froud(num_points))
    allocate(courant(num_points-1))
    allocate(oldQ(num_points, totalChannels))
    allocate(newQ(num_points, totalChannels))
    allocate(oldArea(num_points, totalChannels))
    allocate(newArea(num_points, totalChannels))
    allocate(oldY(num_points, totalChannels))
    allocate(newY(num_points, totalChannels))
    allocate(lateralFlow(num_points,totalChannels))
    allocate(celerity(num_points, totalChannels))
    allocate(diffusivity(num_points, totalChannels))
    allocate(celerity2(num_points))
    allocate(diffusivity2(num_points))
    allocate(eei(num_points))
    allocate(ffi(num_points))
    allocate(exi(num_points))
    allocate(fxi(num_points))
    allocate(qpx(num_points, totalChannels))
    allocate(qcx(num_points))
    allocate(dimensionless_Cr(num_points-1,totalChannels))
    allocate(dimensionless_Fo(num_points-1,totalChannels))
    allocate(dimensionless_Fi(num_points-1,totalChannels))
    allocate(dimensionless_Di(num_points-1,totalChannels))
    allocate(dimensionless_Fc(num_points-1,totalChannels))
    allocate(dimensionless_D(num_points-1,totalChannels))
    allocate(lowerLimitCount(totalChannels))
    allocate(higherLimitCount(totalChannels))
    allocate(currentRoutingNormal(num_points-1,totalChannels))
    allocate(routingNotChanged(num_points-1,totalChannels))
    allocate(elevTable(nel))
    allocate(areaTable(nel))
    allocate(pereTable(nel))
    allocate(rediTable(nel))
    allocate(convTable(nel))
    allocate(topwTable(nel))
    allocate(skkkTable(nel))
    allocate(nwi1Table(nel))
    allocate(dPdATable(nel))
    allocate(ncompElevTable(nel))
    allocate(ncompAreaTable(nel))
    allocate(xsec_tab(11, nel, num_points, totalChannels))
    allocate(rightBank(num_points, totalChannels), leftBank(num_points, totalChannels))
    allocate(skLeft(num_points, totalChannels), skMain(num_points, totalChannels), skRight(num_points, totalChannels))
    allocate(currentSquareDepth(nel))
	allocate(ini_y(nlinks))
	allocate(ini_q(nlinks))
    allocate(notSwitchRouting(nlinks))
    allocate(currentROutingDiffusive(nlinks))

    !dx=dx_ar_g     ! change Nazmul 20210601
    dx = 0.
    minDx = 1e10
    do j = 1,nlinks
        ncomp= frnw_g(j,1)
        do i = 1,ncomp-1
            dx(i,j) = dx_ar_g(i,j)
        end do
        minDx = min(minDx,minval(dx(1:ncomp-1,j)))
    end do

    z=z_ar_g
    ini_y=0.05  !* [meter]
    ini_q=0.1   !*[m^3/sec]
    oldQ = -999; oldY = -999; newQ = -999; newY = -999
    dimensionless_Cr = -999; dimensionless_Fo = -999; dimensionless_Fi = -999
    dimensionless_Di = -999; dimensionless_Fc = -999; dimensionless_D = -999
    ! Reading Q-SK table data data starts
	allocate(noQSKtable(nlinks))
    applyNaturalSection=1
    dt = dtini
    !minDx=minval(dx) ! change Nazmul 20210601

    ! reading bank locations
    if (applyNaturalSection .eq. 1) then
        do j = 1,nlinks
            ncomp= frnw_g(j,1)
            do i=1,ncomp
                leftBank(i,j)= (twcc_ar_g(i,j)-tw_ar_g(i,j))/2.0
                rightBank(i,j)= (twcc_ar_g(i,j)-tw_ar_g(i,j))/2.0 + tw_ar_g(i,j)
            enddo
        enddo
    end if

    do j = 1,nlinks
        ncomp= frnw_g(j,1)
        if (applyNaturalSection .eq. 0) then

        else
            do i=1,ncomp
                skLeft(i,j)= 1.0/manncc_ar_g(i,j)
                skRight(i,j)= 1.0/manncc_ar_g(i,j)
                skMain(i,j)= 1.0/mann_ar_g(i,j)

                call readXsection(i,(1.0/skLeft(i,j)),(1.0/skMain(i,j)),(1.0/skRight(i,j)),&
                    leftBank(i,j), rightBank(i,j),timesDepth,j)

                oldY(i,j) = ini_y(j) + z(i,j)
                oldQ(i,j) = ini_q(j)
            end do
        end if
    end do
    ! reading Q-Strickler's coefficient multiplier table
    do j = 1,nlinks
        ncomp= frnw_g(j,1)
        noQSKtable(i)=0    !* ignore the following table and thus the muliplier is always one.
    end do

    x = 0.0
!**
!* START: v02 routing framework input data
!**
    !* upstream boundary condition data
    allocate( ubcd_g(nts_ub_g, nlinks) )
    do j=1, nlinks
        do n=1,nts_ub_g
            ubcd_g(n,j) =-1.0 !* initialize
        end do
    end do
    !* ubcd_g in [m^3/sec]
    do i1=1, nl_ubcd_ar_g
        read(32,*) frj, n, tarr_ub(n), ubcd_g(n,frj)
    enddo
    !* TW water depth data as downstream boundary data(likely to have missing data)
    allocate(dmyv1(nts_db_g), dmyv2(nts_db_g))
    ndata=0
    do n=1, nts_db_g
        !read(33,*) tarr_db(n), dbcd_g(n)  !* time in [minutes]
        read(33,*) dmy1, dmy2  !* time in [minutes]
        if (dmy1.ge.0.0) then
            ndata=ndata+1
            dmyv1(ndata)= dmy1
            dmyv2(ndata)= dmy2
        end if
    end do
    nts_db_g=ndata
    allocate(tarr_db(nts_db_g), varr_db(nts_db_g))
    allocate(dbcd_g(nts_db_g) )
    tarr_db= dmyv1
    dbcd_g= dmyv2
!**
!* END: v02 routing framework input data
!**

    t=t0*60.0     !! from now on, t is in minute
    !* initial value at initial time for head nodes of head water reach or TW node
    do j = 1, nlinks
        ncomp= frnw_g(j,1)
        if (frnw_g(j,3)==0) then !* frnw_g(j,3) indicates the number of upstream reaches.
        !* head water reach
            do n=1,nts_ub_g
                varr_ub(n)= ubcd_g(n,j)
            enddo
            oldQ(1,j)= intp_y(nts_ub_g, tarr_ub, varr_ub, t) !* tarr_ub in min.
        endif

        if (frnw_g(j,2)<0.0) then
            !* downstream boundary node at TW
            do n=1,nts_db_g
                varr_db(n)= dbcd_g(n) + z(ncomp,j) !* when dbcd is water stage, channel bottom elev is added.
            enddo
            oldY(ncomp,j)= intp_y(nts_db_g, tarr_db, varr_db, t)

            ncompElevTable = xsec_tab(1,:,ncomp,j)
            ncompAreaTable = xsec_tab(2,:,ncomp,j)
            xt=oldY(ncomp,j)

            if (applyNaturalSection .eq. 0) then
                !oldArea(ncomp,j) = ( oldY(ncomp,j) - z(ncomp,j) ) * bo(ncomp,j)
            else
                call r_interpol(ncompElevTable,ncompAreaTable,nel,xt,oldArea(ncomp,j))
                if (oldArea(ncomp,j) .eq. -9999) then
                    print*, 'At j = ',j,', i = ',ncomp, 'time =',t, 'interpolation of oldArea(ncomp,j) was not possible'
                    stop
                end if
            end if
        end if
    end do
    !* correcting the WL initial condition based on the WL boundary
    !* so that the initial WL is higher than or equal to the WL boundary, at j = nlinks, i=ncomp
    do j = 1,nlinks
        ncomp= frnw_g(j,1)
        do i=1,ncomp
            if (oldY(i,j) .lt. oldY(ncomp,nlinks)) oldY(i,j) = oldY(ncomp,nlinks)
        end do
    end do
    ! Applying initial condition of area
    do j=1, nlinks
        ncomp= frnw_g(j,1)
        do i=1,ncomp
            if (applyNaturalSection .eq. 0) then
!                oldArea(i,j) = ( oldY(i,j) - z(i,j) ) * bo(i,j)
            else
                elevTable = xsec_tab(1,1:nel,i,j)
                areaTable = xsec_tab(2,1:nel,i,j)
                call r_interpol(elevTable,areaTable,nel,oldY(i,j),oldArea(i,j))
                if (oldArea(i,j) .eq. -9999) then
                    print*, 'At j = ',j,', i = ',i, 'time =',t, 'interpolation of oldArea(i,j) was not possible'
                    stop
                end if
            end if
        enddo
    end do

    volRemain = -999
    do j=1, nlinks
        ncomp= frnw_g(j,1)
        do i=1,ncomp-1
            volRemain(i,j) = (oldArea(i,j)+oldArea(i+1,j))/2.0*dx(i,j)
        end do
    end do
	 write(101,"(A10, A6, 2A15, 3A15)" ) "time", "nodeID", "reachID", "linkID", "water depth[m]", "Q[cms]", "elevation[m]"
	! Some essential initial parameters for Diffusive Wave
	theta = 1.0
	qpx = 0.
	cfl=0.9
    width = 100. !   initialization
    celerity = 1.0
    maxCelerity = 1.0
    diffusivity = 10.
    maxCelDx = maxCelerity / minDx
    !!! setting initial values of dimensionless parameters
    !dimensionless_Cr, dimensionless_Fo, dimensionless_Fi, dimensionless_Fc, dimensionless_Di, dimensionless_D
    dimensionless_Fi = 10.1
    dimensionless_Fc = 10.1
    dimensionless_D  = 0.1
    currentROutingDiffusive = 1
    ! parameters for diffusive vs partial diffusive
    currentRoutingNormal = 0
    routingNotChanged = 0
    frus2 = 9999.
    notSwitchRouting=0
    minNotSwitchRouting = 10000         ! works between Dynamic and Diffusive switching
    minNotSwitchRouting2 = 000        ! works between full Diffusive and partial Diffusive switching
    timestep = 0


    ! change Nazmul 20210601
    ! Writing initial data
    do j = 1, nlinks
        ncomp= frnw_g(j,1)
        do i=1, ncomp
            write(*,*) t, i, j, oldY(i,j)-z(i,j), oldQ(i,j)
            if (i==1) then
                !* only head node has link ID
                write(101,"(f10.2, i6, 2i15, 3f15.2)" ) t, i, j, pynw(j,2), oldY(i,j)-z(i,j), oldQ(i,j), oldY(i,j)
            else
                idmy1=-100
                write(101,"(f10.2, i6, 2i15, 3f15.2)" ) t, i, j, idmy1, oldY(i,j)-z(i,j), oldQ(i,j), oldY(i,j)
            endif
        enddo
    enddo



    do while ( t .lt. tfin *60.)
        timestep = timestep + 1
        !+-------------------------------------------------------------------------------------
        !+                                      PREDICTOR
        !+
        !+-------------------------------------------------------------------------------------
        do j = 1, nlinks
            ncomp= frnw_g(j,1)
            !+++-- Checking the dtini for possible diffusive wave model and applying it to the model.
            if (j .eq. 1) call calculateDT(t0, t,saveInterval, cfl, tfin, maxCelDx,dtini_given)

            if (frnw_g(j,3)==0) then !* frnw_g(j,3) indicates the number of upstream reaches.
            !* head water reach
                do n=1,nts_ub_g
                    varr_ub(n)= ubcd_g(n,j)
                enddo
                newQ(1,j)= intp_y(nts_ub_g, tarr_ub, varr_ub, t+dtini/60.) !* tarr_ub in min.
            endif
            if (frnw_g(j,2)<0.0) then
                !* downstream boundary node at TW
                do n=1,nts_db_g
                    varr_db(n)= dbcd_g(n) + z(ncomp,j) !* when dbcd is water stage, channel bottom elev is added.
                enddo
                newY(ncomp,j)= intp_y(nts_db_g, tarr_db, varr_db, t+dtini/60.0)

                xt=newY(ncomp,j)

                if (applyNaturalSection .eq. 0) then
    !                 newArea(ncomp,j) = (newY(ncomp,j) - z(ncomp,j)) * bo(ncomp,j)
                else
                    ncompElevTable = xsec_tab(1,:,ncomp,j)
                    ncompAreaTable = xsec_tab(2,:,ncomp,j)
                    call r_interpol(ncompElevTable,ncompAreaTable,nel,xt,newArea(ncomp,j))
                    if (newArea(ncomp,j) .eq. -9999) then
                        print*, 'At j = ',j,', i = ',i, 'time =',t, 'interpolation of newArea(ncomp,j) was not possible'
                        stop
                    end if
                end if
            end if
            !* estimate lateral flow at current time t
            do i=1,ncomp-1
                do n=1,nts_ql_g
                    varr_ql(n)= qlat_g(n,i,j) !* qlat_g(n,i,j) in unit of m2/sec
                enddo
                lateralFlow(i,j)= intp_y(nts_ql_g, tarr_ql, varr_ql, t)
            enddo
            !+++----------------------------------------------------------------
            !+ Hand over water from upstream to downstream properly according
            !+ to the nature of link connections, i.e., serial or branching.
            !+ Refer to p.52,RM1_MESH
            !+++----------------------------------------------------------------
            if (frnw_g(j,3)>0) then !* frnw_g(j,3) indicates the number of upstream reaches.
                !*total water areas at n+1 at the end nodes of upstream links that join link j
                areasum=0.0
                do k=1, frnw_g(j,3)  !* frnw_g(j,3) indicates the number of upstream reaches.
                    linknb= frnw_g(j,3+k)
                    nodenb= frnw_g(linknb,1)
                    areasum=areasum + oldArea(nodenb,linknb) + dap(nodenb,linknb)
                end do
                dqp(1,j)=0.0;
                yav=0.0
                sumOldQ = 0.0
                do k=1, frnw_g(j,3)  !* frnw_g(j,3) indicates the number of upstream reaches.
                    linknb= frnw_g(j,3+k)
                    nodenb= frnw_g(linknb,1)
                    !**dqp(1,j)
                    dqp(1,j)=dqp(1,j)+dqp(nodenb,linknb)    !! Not right! If initial condition is not as sum of Q is conversed, it will be wrong
                    sumOldQ=sumOldQ+oldQ(nodenb,linknb)
                    !**dap(1,j)
                    !*area at the end nod of link k at time n+1
                    areak_ncomp = oldArea(nodenb,linknb) + dap(nodenb,linknb)

                    if (applyNaturalSection .eq. 0) then
    !                    yk_ncomp = areak_ncomp / bo(nodenb,linknb) + z(nodenb,linknb)
                    else
                        elevTable = xsec_tab(1,:,nodenb,linknb)
                        areaTable = xsec_tab(2,:,nodenb,linknb)
                        call r_interpol(areaTable,elevTable,nel,areak_ncomp,yk_ncomp)
                        !* weighted average based on areas at the end nodes of upstream link ks
                    end if
                    yav = yav + (areak_ncomp/areasum)*yk_ncomp
                end do
                dqp(1,j)=dqp(1,j)+sumOldQ-oldQ(1,j) ! Change from DongHa
                newQ(1,j)=oldQ(1,j)+dqp(1,j)        ! If this data is needed in diffusive wave, it takes newQ
                !* Area estimated for time n+1
                if (applyNaturalSection .eq. 0) then
    !                areav = ( yav - z(1,j) ) * bo(1,j)
                else
                    elevTable = xsec_tab(1,:,1,j)
                    areaTable = xsec_tab(2,:,1,j)
                    call r_interpol(elevTable,areaTable,nel,yav,areav)
                end if
                dap(1,j) = areav - oldArea(1,j)
            else        ! There are no links at the upstream of the reach. Example: j = 1, 2
                ! Set upstream discharge
                dqp(1,j) = newQ(1,j) - oldQ(1,j)
                dap(1,j) = 0.0
            end if

            newQ(1,j) = newQ(1,j)+lateralFlow(1,j)*dx(1,j)
            dqp(1,j)  = newQ(1,j) - oldQ(1,j)

           ! checking the value of Fc and Fi in each river reach
            lowerLimitCount = 0; higherLimitCount = 0

            do i=1,ncomp-1
                if ((dimensionless_Fi(i,j) .ge. 5.) .or. (dimensionless_Fc(i,j)  .ge. 5.))then
                    higherLimitCount(j) = higherLimitCount(j) + 1
                elseif ((dimensionless_Fi(i,j) .le. 3.) .or. (dimensionless_Fc(i,j)  .le. 3.))then
                    lowerLimitCount(j) = lowerLimitCount(j) + 1
                end if
            end do
            ! new switching algorithm
            ! for now, auto switching of routing is disabled
            ! manual routing selection:
            ! higherLimitCount(j) = 0; lowerLimitCount(j) = ncomp         ! for dynamic
            ! higherLimitCount(j) = ncomp; lowerLimitCount(j) = ncomp         ! for diffusive
            !! Forcing all computation to diffusive routing
            higherLimitCount = ncomp; lowerLimitCount = ncomp;
            currentROutingDiffusive(j)=1 !* added by DongHa to force diffusive all the time.
            ! running either dynamic or diffusive wave routing at each river reach
            if (higherLimitCount(j) .ge. ncomp/2.) then
                if ( (currentROutingDiffusive(j) .eq. 0) .and. (notSwitchRouting(j) .lt. minNotSwitchRouting)) then
    !                call mesh_dynamic_predictor(dtini_given, t0, t, tfin, saveInterval,j)
                    currentROutingDiffusive(j) = 0
                else
                    call mesh_diffusive_forward(dtini_given, t0, t, tfin, saveInterval,j)
                    if (currentROutingDiffusive(j) .eq. 0) notSwitchRouting(j) = 0
                    currentROutingDiffusive(j) = 1
                end if
            elseif (lowerLimitCount(j) .ge. ncomp/2.) then

                if ( (currentROutingDiffusive(j) .eq. 1) .and. (notSwitchRouting(j) .lt. minNotSwitchRouting)) then
                    call mesh_diffusive_forward(dtini_given, t0, t, tfin, saveInterval,j)
                    currentROutingDiffusive(j) = 1
                else
    !                call mesh_dynamic_predictor(dtini_given, t0, t, tfin, saveInterval,j)
                    if (currentROutingDiffusive(j) .eq. 1) notSwitchRouting(j) = 0
                    currentROutingDiffusive(j) = 0
                end if
            else
                if (currentROutingDiffusive(j) .eq. 1) then
                    call mesh_diffusive_forward(dtini_given, t0, t, tfin, saveInterval,j)
                    currentROutingDiffusive(j) = 1
                else
    !                call mesh_dynamic_predictor(dtini_given, t0, t, tfin, saveInterval,j)
                    currentROutingDiffusive(j) = 0
                end if
            end if

            notSwitchRouting(j) = notSwitchRouting(j) + 1

        end do  ! end off j loop for predictor
        !+-------------------------------------------------------------------------------------
        !+                                      CORRECTOR
        !+
        !+-------------------------------------------------------------------------------------
        do j =  nlinks,1,-1
            ncomp= frnw_g(j,1)
            !+++------------------------------------------------------------+
            !+ Handle downstream boundary condition for a link j that has a link
            !+ immediately downstream.
            !+ **Note. dac and dqc at the last node of the most downstream
            !+ link are computed in sub dsbc during predictor step.
            !+ Refer to p.53-1,RM1_MESH
            !+++------------------------------------------------------------+
            if (frnw_g(j,2).ge.0.0) then !* frnw_g(j,2) gives j of a downstream reach
            !* NOT TW reach
                linknb= frnw_g(j,2)
                newY(ncomp,j)= newY(1,linknb)
                xt = newY(ncomp,j)
                if (applyNaturalSection .eq. 0) then
                    newArea(ncomp,j) = (newY(ncomp,j) - z(ncomp,j)) * bo(ncomp,j)
                else
                    elevTable = xsec_tab(1,:,ncomp,j)
                    areaTable = xsec_tab(2,:,ncomp,j)
                    call r_interpol(elevTable,areaTable,nel,xt,newArea(ncomp,j))
                end if
                dac(ncomp,j)=2*(newArea(ncomp,j)-oldArea(ncomp,j))-dap(ncomp,j)
                areap(ncomp,j) = areap(ncomp,j) - dap(ncomp,j) + dac(ncomp,j)   !! change 20210311 !! new added line
                dqc(ncomp,j)=dqc(1,linknb)*qp(ncomp,j)/qp(1,linknb) ! Changed from what DongHa originally proposed.
                !* p.120,RM3
                qsum= 0.0
                linknb_ds= linknb
                do k=1, frnw_g(linknb_ds,3)
                    linknb_us= frnw_g(linknb_ds,3+k)
                    nodenb= frnw_g(linknb_us,1)
                    qsum= qsum + qp(nodenb,linknb_us)
                end do
                    qnp1_ds= oldQ(1,linknb_ds) +0.5*(dqp(1,linknb_ds)+dqc(1,linknb_ds))
                    !* est. q(n+1, ncomp, link j_i), p120_RM
                    qnp1_us= qnp1_ds*qp(ncomp,j)/qsum
                    dqc(ncomp,j)= 2.0*(qnp1_us - oldQ(ncomp,j)) - dqp(ncomp,j)
                else            !! The river reach is the last reach. i.e. j = 3
                    areap(ncomp,j) = areap(ncomp,j) - dap(ncomp,j) + (newArea(ncomp,j) - oldArea(ncomp,j)) !! change 20210311 !! the calculated areap is now corrected from dac(ncomp)
                    dap(ncomp,j)=newArea(ncomp,j) - oldArea(ncomp,j)    !update from downstream time series
                    dac(ncomp,j)=dap(ncomp,j)
                    dqc(ncomp,j)=dqp(ncomp,j)
                end if

            if (currentROutingDiffusive(j) .eq. 0) then
!                call mesh_dynamic_corrector(dtini_given, t0, t, tfin, saveInterval,j)
            elseif (currentROutingDiffusive(j) .eq. 1) then
                !call mesh_diffusive_backward(dtini_given, t0, t, tfin, saveInterval,j)
                call mesh_diffusive_backward(dtini_given, t0, t, tfin, saveInterval,j,leftBank, rightBank)
            else
                print*, 'Something is wrong in reach ', j
                stop
            end if

            if (j .eq. 1) then
                maxCelDx = 0.
                maxCelerity = 0.
                do i=1,nlinks
                    do kkk = 2, frnw_g(i,1) !nx1(i)
                        maxCelDx = max(maxCelDx,celerity(kkk,i)/dx(kkk-1,i)) ! correction 20210408
                        maxCelerity = max(maxCelerity,celerity(kkk,i))
                    end do
                end do
            end if!
        end do  ! end of j loop

        do j = 1, nlinks
            ncomp= frnw_g(j,1)
            do i=1,ncomp
                froud(i)=abs(newQ(i,j))/sqrt(grav*newArea(i,j)**3.0/bo(i,j))
                if (i .lt. ncomp) then
                    courant(i)=(newQ(i,j)+newQ(i+1,j))/(newArea(i,j)+newArea(i+1,j))*dtini/dx(i,j)
                endif
            enddo
            if (maxCourant .lt. maxval (courant(1:ncomp-1))) then
                maxCourant = maxval (courant(1:ncomp-1))
            endif
        enddo

        do j=1, nlinks
            ncomp= frnw_g(j,1)
            do i=1,ncomp-1
                volRemain(i,j) = (newArea(i,j)+newArea(i+1,j))/2.0*dx(i,j)
            end do
        end do

        t = t + dtini/60.
        ! after a warm up of 24hours, the model will not be forced to run in partial diffusive mode
        if ((t-t0*60.) .ge. 24.*60.) minNotSwitchRouting2 = 100

        do j = 1, nlinks
            ncomp= frnw_g(j,1)
            call calc_dimensionless_numbers(j)
        enddo

        if ( (mod( (t-t0*60.)*60.  ,saveInterval) .le. TOLERANCE) .or. ( t .eq. tfin *60. ) ) then
            do j = 1, nlinks
                ncomp= frnw_g(j,1)
                do i=1, ncomp
                    write(*,*) t, i, j, newY(i,j)-z(i,j), newQ(i,j)
                    if (i==1) then
                        !* only head node has link ID
                        write(101,"(f10.2, i6, 2i15, 3f15.2)" ) t, i, j, pynw(j,2), newY(i,j)-z(i,j), newQ(i,j), newY(i,j)
                    else
                        idmy1=-100
                        write(101,"(f10.2, i6, 2i15, 3f15.2)" ) t, i, j, idmy1, newY(i,j)-z(i,j), newQ(i,j), newY(i,j)
                    endif
                enddo
            enddo
        end if
        ! update of Y, Q and Area vectors
        oldY   = newY
        newY=-999
        oldQ   = newQ
        newQ=-999
        oldArea= newArea
        newArea=-999
        pere=-999
    enddo  ! end of time loop

    call cpu_time( t2 )
    print '("Total CPU Time = ",f10.3," seconds.")',t2 - t1

end program mesh

real function r_interpol_time(x,y,jj,xt)
    implicit none
    integer, intent(in) :: jj
    real, intent(in) :: x(jj), y(jj)
    real, intent(in) :: xt
    real :: yt
    integer :: i, j

    if (xt.le.maxval(x) .and. xt.ge.minval(x)) then
        do j=1,jj-1
            if((x(j)-xt)*(x(j+1)-xt).le.0)then

                yt=(xt-x(j))/(x(j+1)-x(j))*(y(j+1)-y(j))+y(j)

                EXIT
            endif
        end do
    else
        print*, xt, ' is not within the limit'
        print*, 'maxval(x)= ', maxval(x), 'and minval(x)=', minval(x),'so',  xt, ' is not within the limit'
        print*, 'jj', jj
        print*, 'x', (x(i), i=1, jj)
        print*, 'y', (y(i), i=1, jj)
        stop
    end if
    r_interpol_time = yt
    return
end function

real function r_interpo_nn(x,y,jj,xt)
    implicit none
    integer, intent(in) :: jj
    real, intent(in) :: xt, x(jj), y(jj)
    real :: yt
    real :: j
    ! nn means nearest neighbour
    if (xt.le. x(1)) then
        yt=y(1)
    elseif (xt.ge. x(jj)) then
        yt=y(jj)
    else
        do j=1,jj-1
            if((x(j)-xt)*(x(j+1)-xt).le.0)then

                yt=(xt-x(j))/(x(j+1)-x(j))*(y(j+1)-y(j))+y(j)

                EXIT
            endif
        end do
    end if
    r_interpo_nn = yt
    return
end function
