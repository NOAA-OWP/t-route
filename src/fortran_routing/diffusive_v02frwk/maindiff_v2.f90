module mdiffv2    
    !*-------------------------------------------------------------------------------------------------
    !*       A diffusive model developed by Tulane University (Prof.Ehab Meselhe and 
    !*		Dr.Md Nazmul Azim Beg) and integrated into a Python routing framework By Dong Ha Kim at
    !*          the National Water Center. Basically, the partial differential equation of diffusive 
    !*          wave is numerically solved using Crank-Nicoloson and Hermite Interpolation techniques.
    !*          Water depth computation produces either normal depth or diffusive depth, the selection
    !*          of which is determined mainly by dimensionless diffusion coefficient. 
    !*-------------------------------------------------------------------------------------------------

    implicit none
!    !* symbolic names for kind types of 4-, 2-, and 1-byte integers:
!    integer, parameter :: i4b = selected_int_kind(9)
!    integer, parameter :: i2b = selected_int_kind(4)
!    integer, parameter :: i1b = selected_int_kind(2)
!    !* symbolic names for kind types of single- and double-precision reals:
!    integer, parameter :: sp = kind(1.0)
!    integer, parameter :: dp = kind(1.0d0)

    double precision, parameter :: grav = 9.81
    double precision, parameter :: TOLERANCE = 1e-8
    integer :: nlinks, mxncomp, maxTableLength, nel
    double precision :: dtini, dxini, cfl, minDx, maxCelerity,  theta
    double precision :: frus2, minNotSwitchRouting, minNotSwitchRouting2

    double precision, dimension(:), allocatable :: area, depth, co, froud, courant
    double precision, dimension(:,:), allocatable :: bo, dx
	!**arrays for branching channel application
    double precision, dimension(:,:), allocatable :: areap, qp, z, sk
    double precision, dimension(:,:), allocatable :: dqp, dap, dqc, dac
    double precision, dimension(:,:), allocatable :: celerity, diffusivity, qpx
    double precision, dimension(:), allocatable :: eei, ffi, exi, fxi, qcx, diffusivity2, celerity2
    ! change for unsteady flow
    double precision, dimension(:,:), allocatable :: pere, oldQ, newQ, oldArea, newArea, oldY, newY
    double precision, dimension(:,:), allocatable :: lateralFlow
    double precision, dimension(:,:), allocatable :: dimensionless_Cr, dimensionless_Fo, dimensionless_Fi
    double precision, dimension(:,:), allocatable :: dimensionless_Di, dimensionless_Fc, dimensionless_D
    double precision, dimension(:), allocatable :: ini_y, ini_q
    double precision, dimension(:), allocatable :: lowerLimitCount, higherLimitCount
    double precision, dimension(:,:), allocatable :: volRemain
    integer, dimension(:), allocatable :: currentROutingDiffusive, notSwitchRouting
    integer, dimension(:,:), allocatable :: currentRoutingNormal, routingNotChanged

    double precision, dimension(:), allocatable :: elevTable, areaTable, skkkTable
    double precision, dimension(:), allocatable :: pereTable, rediTable
    double precision, dimension(:), allocatable :: convTable, topwTable
    double precision, dimension(:), allocatable :: nwi1Table, dPdATable
    double precision, dimension(:), allocatable :: ncompElevTable, ncompAreaTable
    double precision, dimension(:,:,:,:), allocatable :: xsec_tab
    double precision, dimension(:), allocatable :: currentSquareDepth

    integer, dimension(:,:), allocatable :: frnw_g
    double precision :: z_g, bo_g, traps_g, tw_g, twcc_g, so_g, mann_g, manncc_g
    integer :: applyNaturalSection
    integer :: nel_g

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

        double precision, dimension(nrch_g, frnw_col), intent(in) :: dfrnw_g !* set frnw_col=10
        double precision, dimension(nts_ql_g, mxncomp_g, nrch_g), intent(in) :: qlat_g
        double precision, dimension(nts_ub_g, nrch_g), intent(in) :: ubcd_g
        double precision, dimension(nts_db_g,6), intent(in) :: dbcd_g

        double precision, intent(in) :: cfl_g, theta_g, so_llm_g
        integer, intent(in) :: tzeq_flag_g !* 0 for lookup tabale; 1 for using procedures to compute trapz.ch.geo.
        integer, intent(in) :: y_opt_g  !* 1 for normal depth(kinematic); 2 for dept of diffusive wave.

        double precision, dimension(mxncomp_g, nrch_g), intent(inout) :: so_ar_g
        double precision, dimension(ntss_ev_g, mxncomp_g, nrch_g), intent(out) :: q_ev_g, elv_ev_g
        integer :: ncomp

        integer :: i, j, k, ppn, qqn, n, ntim, igate, pp, boundaryFileMaxEntry, saveFrequency
        integer :: linknb_ds, linknb_us
        double precision :: qnp1_ds, qnp1_us, qsum, y_ds
        double precision :: cour, da, dq, x, saveInterval, width
        double precision :: qn, xt, maxCourant, dtini_given, nodenb, linknb
        double precision :: frds, areasum, yk_ncomp, yav, areak_ncomp, areav, sumOldQ, currentQ, area_ds
        double precision :: arean, areac, hyrdn, hyrdc, perimn, perimc, qcrit, s0ds, timesDepth
        doubleprecision :: latFlowValue, latFlowValue2
        double precision :: t, r_interpol_time, tfin, t1, t2, t0 !t0 start time
        integer :: tableLength, timestep, kkk
        double precision :: area_0, width_0, errorY, hydR_0, q_sk_multi, sumCelerity
        double precision :: r_interpo_nn
        double precision :: maxCelDx
        double precision, dimension(:), allocatable :: dmyv, dmyv1, dmyv2
        double precision, dimension(:), allocatable ::tarr_ql, varr_ql, tarr_ub, varr_ub, tarr_db, varr_db
        integer :: frj, iseg, i1, ts_ev
        integer :: num_points, totalChannels
        double precision, dimension(:,:), allocatable :: leftBank, rightBank
        double precision, dimension(:,:), allocatable :: skLeft, skMain, skRight
        double precision :: dmy1, dmy2
        integer :: ndata, idmy1, nts_db_g2

        !open(unit=101, file="./output/simulated discharge depth elev.txt")

        nlinks=nrch_g
        allocate(frnw_g(nlinks,frnw_col))
        frnw_g=dfrnw_g
        mxncomp= mxncomp_g
        dtini= dtini_g
        dtini_given= dtini
        t0=t0_g
        tfin= tfin_g
        ntim = floor( (tfin - t0) / dtini * 3600)
        timesDepth= 4.0 !* water depth multiplier used in readXsection
        nel= nel_g
        saveInterval= saveInterval_ev_g
        saveFrequency = saveInterval / dtini_given
        num_points= mxncomp
        totalChannels= nlinks
        allocate(area(num_points))
        ! change for unsteady flow
        allocate(bo(num_points,totalChannels))
        allocate(pere(num_points,totalChannels))
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
        !allocate(dt(num_points))
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
        allocate(tarr_ql(nts_ql_g), varr_ql(nts_ql_g))
        allocate(tarr_ub(nts_ub_g), varr_ub(nts_ub_g))

        !* for dbcd_g, find the actual length of array without missing data
        nts_db_g2=0
        for i=1, nts_db_g
            if (dbcd_g(i).ge.0.0) then
                nts_db_g2= nts_db_g2+1
            end if
        endif
        allocate(tarr_db(nts_db_g2), varr_db(nts_db_g2))

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
        !* Reading Q-SK table data data starts
        applyNaturalSection=1
        !* reading bank locations
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
                                        leftBank(i,j), rightBank(i,j),timesDepth, j,&
                                        z_ar_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g)

                    oldY(i,j) = ini_y(j) + z(i,j)
                    oldQ(i,j) = ini_q(j)
                end do
            end if
        end do
        ! reading Q-Strickler's coefficient multiplier table
        do j = 1,nlinks
            ncomp= frnw_g(j,1)
            !noQSKtable(i)=0    !* ignore the following table and thus the muliplier is always one.
        end do
        x = 0.0

        !* time step series for lateral flow
        do n=1, nts_ql_g
            tarr_ql(n)= t0_g*60.0 + dt_ql_g*real(n-1,KIND(dt_ql_g))/60.0 !* [min]
        end do
        !* time step series for upstream boundary data
        do n=1, nts_ub_g
            tarr_ub(n)= t0_g*60.0 + dt_ub_g*real(n-1,KIND(dt_ub_g))/60.0 !* [min]
        enddo
        !* time step series for downstream boundary data
        do n=1, nts_db_g2
            tarr_db(n)= dbcd_g(n,1)         !* [min]
            !tarr_db(n)= dbcd_g(n,1) t0_g*60.0 + dt_db_g*real(n-1,KIND(dt_db_g))/60.0
        enddo

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
                do n=1,nts_db_g2
                    varr_db(n)= dbcd_g(n,2) + z(ncomp,j) !* when dbcd is water stage, channel bottom elev is added.
                enddo
                oldY(ncomp,j)= intp_y(nts_db_g2, tarr_db, varr_db, t)

                ncompElevTable = xsec_tab(1,:,ncomp,j)
                ncompAreaTable = xsec_tab(2,:,ncomp,j)
                xt=oldY(ncomp,j)

                if (applyNaturalSection .eq. 0) then
                    oldArea(ncomp,j) = ( oldY(ncomp,j) - z(ncomp,j) ) * bo(ncomp,j)
                else
                    call r_interpol(ncompElevTable,ncompAreaTable,nel,xt,oldArea(ncomp,j))
                    if (oldArea(ncomp,j) .eq. -9999) then
                        !print*, 'At j = ',j,', i = ',ncomp, 'time =',t, 'interpolation of oldArea(ncomp,j) was not possible'
                        !stop
                    end if
                end if
            end if
        enddo
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
                    oldArea(i,j) = ( oldY(i,j) - z(i,j) ) * bo(i,j)
                else
                    elevTable = xsec_tab(1,1:nel,i,j)
                    areaTable = xsec_tab(2,1:nel,i,j)
                    call r_interpol(elevTable,areaTable,nel,oldY(i,j),oldArea(i,j))
                    if (oldArea(i,j) .eq. -9999) then
                        !print*, 'At j = ',j,', i = ',i, 'time =',t, 'interpolation of oldArea(i,j) was not possible'
                        !stop
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
         !write(101,"(A10, A6, 2A15, 3A15)" ) "time", "nodeID", "reachID", "linkID", "water depth[m]", "Q[cms]", "elevation[m]"
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
        ts_ev=1 !* time step for outputting q and elv at evaluation time

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
                    do n=1,nts_db_g2
                        varr_db(n)= dbcd_g(n) + z(ncomp,j) !* when dbcd is water stage, channel bottom elev is added.
                    enddo
                    newY(ncomp,j)= intp_y(nts_db_g2, tarr_db, varr_db, t+dtini/60.0)

                    xt=newY(ncomp,j)

                    if (applyNaturalSection .eq. 0) then
                        newArea(ncomp,j) = (newY(ncomp,j) - z(ncomp,j)) * bo(ncomp,j)
                    else
                        ncompElevTable = xsec_tab(1,:,ncomp,j)
                        ncompAreaTable = xsec_tab(2,:,ncomp,j)
                        call r_interpol(ncompElevTable,ncompAreaTable,nel,xt,newArea(ncomp,j))
                        if (newArea(ncomp,j) .eq. -9999) then
!                            print*, 'At j = ',j,', i = ',i, 'time =',t, 'interpolation of newArea(ncomp,j) was not possible'
!                            stop
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
                    !* total water areas at n+1 at the end nodes of upstream links that join link j
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
                            yk_ncomp = areak_ncomp / bo(nodenb,linknb) + z(nodenb,linknb)
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
                        areav = ( yav - z(1,j) ) * bo(1,j)
                    else
                        elevTable = xsec_tab(1,:,1,j)
                        areaTable = xsec_tab(2,:,1,j)
                        call r_interpol(elevTable,areaTable,nel,yav,areav)
                    end if
                    dap(1,j) = areav - oldArea(1,j)
                else        ! There are no links at the upstream of the reach. Example: j = 1, 2
                    !* Set upstream discharge
                    dqp(1,j) = newQ(1,j) - oldQ(1,j)
                    dap(1,j) = 0.0
                end if

                newQ(1,j) = newQ(1,j)+lateralFlow(1,j)*dx(1,j)
                dqp(1,j)  = newQ(1,j) - oldQ(1,j)

               !* checking the value of Fc and Fi in each river reach
                lowerLimitCount = 0; higherLimitCount = 0

                do i=1,ncomp-1
                    if ((dimensionless_Fi(i,j) .ge. 5.) .or. (dimensionless_Fc(i,j)  .ge. 5.))then
                        higherLimitCount(j) = higherLimitCount(j) + 1
                    elseif ((dimensionless_Fi(i,j) .le. 3.) .or. (dimensionless_Fc(i,j)  .le. 3.))then
                        lowerLimitCount(j) = lowerLimitCount(j) + 1
                    end if
                end do
                !** new switching algorithm
                !* for now, auto switching of routing is disabled
                !* manual routing selection:
                !* For dynamic, higherLimitCount(j) = 0; lowerLimitCount(j) = ncomp
                !* For diffusive, higherLimitCount(j) = ncomp; lowerLimitCount(j) = ncomp
                !* Forcing all computation to diffusive routing
                higherLimitCount = ncomp; lowerLimitCount = ncomp;
                currentROutingDiffusive(j)=1 !* added by DongHa to force diffusive all the time.
                !* running either dynamic or diffusive wave routing at each river reach
                if (higherLimitCount(j) .ge. ncomp/2.) then
                    if ( (currentROutingDiffusive(j) .eq. 0) .and. (notSwitchRouting(j) .lt. minNotSwitchRouting)) then
                        !call mesh_dynamic_predictor(dtini_given, t0, t, tfin, saveInterval,j)
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
                        !call mesh_dynamic_predictor(dtini_given, t0, t, tfin, saveInterval,j)
                        if (currentROutingDiffusive(j) .eq. 1) notSwitchRouting(j) = 0
                        currentROutingDiffusive(j) = 0
                    end if
                else
                    if (currentROutingDiffusive(j) .eq. 1) then
                        call mesh_diffusive_forward(dtini_given, t0, t, tfin, saveInterval,j)
                        currentROutingDiffusive(j) = 1
                    else
                        !call mesh_dynamic_predictor(dtini_given, t0, t, tfin, saveInterval,j)
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
                else
                !* TW reach
                    areap(ncomp,j) = areap(ncomp,j) - dap(ncomp,j) + (newArea(ncomp,j) - oldArea(ncomp,j)) !! change 20210311 !! the calculated areap is now corrected from dac(ncomp)
                    dap(ncomp,j)=newArea(ncomp,j) - oldArea(ncomp,j)    !update from downstream time series
                    dac(ncomp,j)=dap(ncomp,j)
                    dqc(ncomp,j)=dqp(ncomp,j)
                end if

                if (currentROutingDiffusive(j) .eq. 0) then
                    !call mesh_dynamic_corrector(dtini_given, t0, t, tfin, saveInterval,j)
                elseif (currentROutingDiffusive(j) .eq. 1) then
                    call mesh_diffusive_backward(dtini_given, t0, t, tfin, saveInterval,j,leftBank, rightBank)
                else
                    !print*, 'Something is wrong in reach ', j
                    !stop
                end if

                if (j .eq. 1) then
                    maxCelDx = 0.
                    maxCelerity = 0.
                    do i=1,nlinks
                        do kkk = 2, frnw_g(i,1)
                            maxCelDx = max(maxCelDx,celerity(kkk,i)/dx(kkk-1,i)) ! correction 20210408
                            maxCelerity = max(maxCelerity,celerity(kkk,i))
                        end do
                    end do
                endif
            enddo  ! end of j loop

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
            !* after a warm up of 24hours, the model will not be forced to run in partial diffusive mode
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
                        q_ev_g(ts_ev, i, j)= newQ(i,j)
                        elv_ev_g(ts_ev, i, j)= newY(i,j)
                        !idmy1=-100
                        !write(101,"(f10.2, i6, 2i15, 3f15.2)" ) t, i, j, idmy1, newY(i,j)-z(i,j), newQ(i,j), newY(i,j)
                    enddo
                enddo
                ts_ev=ts_ev+1
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
    endsubroutine diffnw
    !*--------------------------------------------
    !          Interpolate in time
    !
    !*--------------------------------------------
    double precision function r_interpol_time(x,y,jj,xt)
        implicit none
        integer, intent(in) :: jj
        doubleprecision, intent(in) :: x(jj), y(jj)
        doubleprecision, intent(in) :: xt
        doubleprecision :: yt
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
    endfunction
    !*-----------------------------------------------
    !          Interpolate by nearest neighbor
    !
    !*-----------------------------------------------
    double precision function r_interpo_nn(x,y,jj,xt)
        implicit none
        integer, intent(in) :: jj
        double precision, intent(in) :: xt, x(jj), y(jj)
        double precision :: yt
        double precision :: j
        ! nn means nearest neighbor
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
    !*----------------------------------------------------------------------------------
    !               Adjust simulation time step according to Curant Condition
    !
    !*----------------------------------------------------------------------------------
    subroutine calculateDT(initialTime, time, saveInterval, maxAllowCourantNo, tfin, max_C_dx, given_dt)
        implicit none

        double precision, intent(in) :: initialTime, time, saveInterval, tfin, given_dt
        double precision, intent(in) :: maxAllowCourantNo, max_C_dx
        integer          :: a, b

        !! initialTime is in hours
        !! tfin is in hours
        !! time is in minutes
        !! dtini is in seconds
        !! saveInterval is in seconds
        dtini = maxAllowCourantNo/max_C_dx

        a = floor( (time-initialTime*60.) /( saveInterval/60. ))            ! units:: time : min;  ! initialTime : hour ! saveInterval : sec
        b = floor( (time-initialTime*60.) +dtini/60.)/( saveInterval/60. )

        if (b .gt. a) then
            dtini = (a+1) * ( saveInterval ) - (time-initialTime*60.)*60.
        end if

        if ( time+dtini/60. .gt. tfin*60. ) dtini =  (tfin*60.-time)*60.
    end subroutine
    !**--------------------------------------------------------------------------------
    !*      Compute dimensionless parameters mainly for deciding which depth
    !*      computation schemes to use between normal depth and diffusive depth.
    !
    !**--------------------------------------------------------------------------------
    subroutine calc_dimensionless_numbers(j)
        implicit none

        integer, intent(in) :: j
        integer :: i, ncomp
        double precision :: wl_us, depth_us, q_us, v_us, pere_us, r_us, sk_us, ch_us
        double precision :: wl_ds, depth_ds, q_ds, v_ds, pere_ds, r_ds, sk_ds, ch_ds
        double precision :: ch_star_avg, channel_length, avg_celerity, avg_velocity, avg_depth
        double precision :: maxValue, dimlessWaveLength
        maxValue = 1e7

        ncomp= frnw_g(j,1)
        do i=1, ncomp-1

            avg_celerity = (celerity(i,j) + celerity(i+1,j)) / 2.0    ! 'celerity2' is calculated celerity. 'celerity' is spatially averaged

            wl_us = newY(i,j)
            ! bo(i) and pere(i) has data for the latest river reach only
            depth_us = newArea(i,j) / bo(i,j)
            q_us = newQ(i,j)
            v_us = abs( newQ(i,j) / newArea(i,j) )
            pere_us = pere(i,j)
            r_us = newArea(i,j) / pere(i,j)
            sk_us = sk(i,j)
            ch_us = sk(i,j) * r_us ** (1./6.)


            wl_ds = newY(i+1,j)
            depth_ds = newArea(i+1,j) / bo(i+1,j)
            q_ds = newQ(i+1,j)
            v_ds = abs( newQ(i+1,j) / newArea(i+1,j) )
            pere_ds = pere(i+1,j)
            r_ds = newArea(i+1,j) / pere(i+1,j)
            sk_ds = sk(i+1,j)
            ch_ds = sk(i+1,j) * r_ds ** (1./6.)


            ch_star_avg = ((ch_us + ch_ds) / 2.)  / sqrt( grav ) !! CORRECTED
            channel_length = dx(i,j)

            dimlessWaveLength= 4000. !! new

            avg_velocity = (v_us + v_ds) / 2.
            avg_depth = (depth_us + depth_ds) / 2.

            dimensionless_Cr(i,j) = abs(avg_velocity / avg_celerity)
            if (dimensionless_Cr(i,j) .gt. maxValue) dimensionless_Cr(i,j) = maxValue

            dimensionless_Fo(i,j) = avg_velocity / sqrt(grav * avg_depth)
            if (dimensionless_Fo(i,j) .gt. maxValue) dimensionless_Fo(i,j) = maxValue

            dimensionless_Fi(i,j) = 2*dimensionless_Cr(i,j) / (ch_star_avg ** 2.) * dimlessWaveLength !(channel_length / avg_depth) !! CORRECTED
            if (dimensionless_Fi(i,j) .gt. maxValue) dimensionless_Fi(i,j) = maxValue

            dimensionless_Fc(i,j) = dimensionless_Cr(i,j) * dimensionless_Fi(i,j)
            if (dimensionless_Fc(i,j) .gt. maxValue) dimensionless_Fc(i,j) = maxValue

            dimensionless_Di(i,j) = (dimensionless_Cr(i,j) / dimensionless_Fo(i,j)) ** 2. !! CORRECTED
            if (dimensionless_Di(i,j) .gt. maxValue) dimensionless_Di(i,j) = maxValue

            dimensionless_D(i,j)  = dimensionless_Di(i,j) / dimensionless_Fc(i,j)
            if (dimensionless_D(i,j) .gt. maxValue) dimensionless_D(i,j) = maxValue
        end do
    end subroutine calc_dimensionless_numbers
    !**--------------------------------------------------------------------------------
    !*      Compute discharge using diffusive equation that is numerically solved by
    !*      Crank-Nicolson + Hermite Interpolation method
    !
    !**--------------------------------------------------------------------------------
    subroutine mesh_diffusive_forward(dtini_given, t0, t, tfin, saveInterval,j)
        implicit none
        integer, intent(in) :: j
        double precision, intent(in) :: dtini_given, t0, t, tfin, saveInterval
        double precision :: a1, a2, a3, a4, b1, b2, b3, b4, dd1, dd2, dd3, dd4, h1, h2, h3, h4, xt, allqlat
        double precision :: qy, qxy, qxxy, qxxxy, ppi, qqi, rri, ssi, sxi, mannings, Sb, width
        double precision :: cour, cour2, q_sk_multi, sfi, temp, alpha !r_interpol_time, r_interpo_nn,
        double precision :: y_norm_ds, y_crit_ds, S_ncomp, frds, area_0, width_0, hydR_0, errorY, currentQ
        integer :: tableLength, ll, ncomp
        integer :: i, pp
        double precision :: eei_ghost, ffi_ghost, exi_ghost, fxi_ghost, qp_ghost, qpx_ghost

        !* change 20210228: All qlat to a river reach is applied to the u/s boundary
        !* Note: lateralFlow(1,j) is already added to the boundary
        ncomp= frnw_g(j,1)
        allqlat = sum(lateralFlow(2:ncomp-1,j) * dx(2:ncomp-1,j)) ! change Nazmul 20210601
        eei = -999.
        ffi = -999. !! What will be this value?
        exi = -999.
        fxi = -999.
        !* steps for advection equation
        eei(1) = 1.0
        ffi(1) = 0. !! What will be this value?
        exi(1) = 0.
        fxi(1) = 0.

        ncomp= frnw_g(j,1)
        do i = 2,ncomp
            !!!------ Calculation a1...a4, up to h4...
            cour = dtini / dx(i-1,j)
            cour2= abs( celerity(i,j) ) * cour

            a1 = 3.0 * cour2 ** 2.0 - 2.0 * cour2 ** 3.0
            a2 = 1 - a1
            a3 = ( cour2 ** 2.0 - cour2 ** 3.0 ) * dx(i-1,j)
            a4 = ( -1.0 * cour2 + 2.0 * cour2 ** 2.0 - cour2 ** 3.0 ) * dx(i-1,j)

            b1 = ( 6.0 * cour2 - 6.0 * cour2 ** 2.0 ) / ( -1.0 * dx(i-1,j) )
            b2 = - b1
            b3 = ( 2.0 * cour2 - 3.0 * cour2 ** 2.0 ) * ( -1.0 )
            b4 = ( -1.0 + 4.0 * cour2 - 3.0 * cour2 ** 2.0 ) * ( -1.0 )

            dd1 = ( 6.0 - 12.0 * cour2 ) / ( dx(i-1,j) ** 2.0 )
            dd2 = - dd1
            dd3 = ( 2.0 - 6.0 * cour2 ) / dx(i-1,j)
            dd4 = ( 4.0 - 6.0 * cour2 ) / dx(i-1,j)

            h1 = 12.0 / ( dx(i-1,j) ** 3.0 )
            h2 = - h1
            h3 = 6.0 / ( dx(i-1,j) ** 2.0 )
            h4 = h3

            if (i .eq. ncomp) then
                alpha = 1.0
            else
                alpha = dx(i,j) / dx(i-1,j)
            end if

            qy   = a1 * oldQ(i-1,j) + a2 * oldQ(i,j) + a3 * qpx(i-1,j) + a4 * qpx(i,j)
            qxy  = b1 * oldQ(i-1,j) + b2 * oldQ(i,j) + b3 * qpx(i-1,j) + b4 * qpx(i,j)
            qxxy = dd1* oldQ(i-1,j) + dd2* oldQ(i,j) + dd3* qpx(i-1,j) + dd4* qpx(i,j)
            qxxxy= h1 * oldQ(i-1,j) + h2 * oldQ(i,j) + h3 * qpx(i-1,j) + h4 * qpx(i,j)

            ppi = - theta * diffusivity(i,j) * dtini / ( dx(i-1,j) ** 2.0 ) * 2.0 / (alpha*(alpha + 1.0)) * alpha
            qqi = 1.0 - ppi * (alpha + 1.0) / alpha
            rri = ppi / alpha

            ssi = qy  + dtini * diffusivity(i,j) * ( 1.0 - theta ) * qxxy
            sxi = qxy + dtini * diffusivity(i,j) * ( 1.0 - theta ) * qxxxy

            eei(i) = -1.0 * rri / ( ppi * eei(i-1) + qqi )                     !! copied from split operator method
            ffi(i) = ( ssi - ppi * ffi(i-1) ) / ( ppi * eei(i-1) + qqi )       !! copied from split operator method
            exi(i) = -1.0 * rri / ( ppi * exi(i-1) + qqi )
            fxi(i) = ( sxi - ppi * fxi(i-1) ) / ( ppi * exi(i-1) + qqi )
        end do
        !!! Ghost point calculation start
        cour = dtini / dx(ncomp-1,j)
        cour2= abs( celerity(ncomp-1,j) ) * cour

        a1 = 3.0 * cour2 ** 2.0 - 2.0 * cour2 ** 3.0
        a2 = 1 - a1
        a3 = ( cour2 ** 2.0 - cour2 ** 3.0 ) * dx(ncomp-1,j)
        a4 = ( -1.0 * cour2 + 2.0 * cour2 ** 2.0 - cour2 ** 3.0 ) * dx(ncomp-1,j)

        b1 = ( 6.0 * cour2 - 6.0 * cour2 ** 2.0 ) / ( -1.0 * dx(ncomp-1,j) )
        b2 = - b1
        b3 = ( 2.0 * cour2 - 3.0 * cour2 ** 2.0 ) * ( -1.0 )
        b4 = ( -1.0 + 4.0 * cour2 - 3.0 * cour2 ** 2.0 ) * ( -1.0 )

        dd1 = ( 6.0 - 12.0 * cour2 ) / ( dx(ncomp-1,j) ** 2.0 )
        dd2 = - dd1
        dd3 = ( 2.0 - 6.0 * cour2 ) / dx(ncomp-1,j)
        dd4 = ( 4.0 - 6.0 * cour2 ) / dx(ncomp-1,j)

        h1 = 12.0 / ( dx(ncomp-1,j) ** 3.0 )
        h2 = - h1
        h3 = 6.0 / ( dx(ncomp-1,j) ** 2.0 )
        h4 = h3

        alpha = 1.0

        qy   = a1 * oldQ(ncomp,j) + a2 * oldQ(ncomp-1,j) + a3 * qpx(ncomp,j) + a4 * qpx(ncomp-1,j)
        qxy  = b1 * oldQ(ncomp,j) + b2 * oldQ(ncomp-1,j) + b3 * qpx(ncomp,j) + b4 * qpx(ncomp-1,j)
        qxxy = dd1* oldQ(ncomp,j) + dd2* oldQ(ncomp-1,j) + dd3* qpx(ncomp,j) + dd4* qpx(ncomp-1,j)
        qxxxy= h1 * oldQ(ncomp,j) + h2 * oldQ(ncomp-1,j) + h3 * qpx(ncomp,j) + h4 * qpx(ncomp-1,j)

        ppi = - theta * diffusivity(ncomp,j) * dtini / ( dx(ncomp-1,j) ** 2.0 ) * 2.0 / (alpha*(alpha + 1.0)) * alpha
        qqi = 1.0 - ppi * (alpha + 1.0) / alpha
        rri = ppi / alpha

        ssi = qy  + dtini * diffusivity(ncomp-1,j) * ( 1.0 - theta ) * qxxy
        sxi = qxy + dtini * diffusivity(ncomp-1,j) * ( 1.0 - theta ) * qxxxy

        eei_ghost = -1.0 * rri / ( ppi * eei(ncomp) + qqi )                     !! copied from split operator method
        ffi_ghost = ( ssi - ppi * ffi(ncomp) ) / ( ppi * eei(ncomp) + qqi )       !! copied from split operator method

        exi_ghost = -1.0 * rri / ( ppi * exi(ncomp) + qqi )
        fxi_ghost = ( sxi - ppi * fxi(ncomp) ) / ( ppi * exi(ncomp) + qqi )

        !!! Ghost point calculation end
        qp_ghost = oldQ(ncomp-1,j)
        qpx_ghost= 0.

        qp(ncomp,j) = eei(ncomp) * qp_ghost + ffi(ncomp)
        qpx(ncomp,j)= exi(ncomp) *qpx_ghost + fxi(ncomp)

        do i = ncomp-1,1,-1
            qp(i,j) = eei(i) * qp(i+1,j) + ffi(i)
            qpx(i,j)= exi(i) *qpx(i+1,j) + fxi(i)
        end do

        qp(1,j) = newQ(1,j)     ! change Nazmul 20210601

        ! change 20210228: All qlat to a river reach is applied to the u/s boundary
        qp(1,j) = qp(1,j) + allqlat

        do i=1,ncomp
            if (abs(qp(i,j)) .lt. 0.02831) then
                qp(i,j) = 0.02831
            end if
        end do

        newQ(1:ncomp,j) = qp(1:ncomp,j)
        dqp(1:ncomp,j) = newQ(1:ncomp,j)-oldQ(1:ncomp,j)
        dqc(1:ncomp,j) = dqp(1:ncomp,j)
        dap(1:ncomp,j) = 0.

    end subroutine mesh_diffusive_forward
    !**--------------------------------------------------------------------------------
    !*      Compute water depth using diffusive momentum equation or normal depth
    !*      with computed Q from the forward step.
    !
    !**--------------------------------------------------------------------------------
    subroutine mesh_diffusive_backward(dtini_given, t0, t, tfin, saveInterval,j, leftBank, rightBank)

        implicit none

        integer, intent(in) :: j
        double precision, intent(in) :: dtini_given, t0, t, tfin, saveInterval
        double precision, dimension(mxncomp, nlinks) :: leftBank, rightBank
        double precision :: a1, a2, a3, a4, b1, b2, b3, b4, dd1, dd2, dd3, dd4, h1, h2, h3, h4, xt
        double precision :: qy, qxy, qxxy, qxxxy, ppi, qqi, rri, ssi, sxi, mannings, Sb, width, slope
        double precision :: cour, cour2, q_sk_multi, sfi, temp, dkdh !r_interpol_time, r_interpo_nn,
        double precision :: D_lim1, D_lim2, y_norm, y_crit, area_n, area_c, chnWidth, vel
        double precision :: y_norm_ds, y_crit_ds, S_ncomp, frds, area_0, width_0, hydR_0, errorY, currentQ, stg1, stg2
        integer :: tableLength, jj, newMassBalance, iii
        double precision :: elevTable_1(nel),areaTable_1(nel),rediTable_1(nel),convTable_1(nel),topwTable_1(nel)
        double precision :: skkkTable_1(nel), dKdATable_1(nel), pereTable_1(nel),currentSquaredDepth_1(nel)                                                                       ! change Nazmul 20210601
        double precision :: depthYi,tempDepthi_1,tempCo_1,temp_q_sk_multi_1,tempY_1,tempArea_1,tempRadi_1,tempbo_1,ffy
        double precision :: ffy_1, ffy_2, ffy_3, tempCo_2, tempCo_3, tempsfi_2, tempsfi_3
        double precision :: ffprime,tempDepthi_1_new,tempsfi_1,toll, dkda
        doubleprecision :: tempPere_1, tempsk_1, tempY_2, tempY_3, tempdKdA_1               ! change Nazmul 20210601
        integer :: depthCalOk(mxncomp), newtonRaphson
        integer :: i, pp, ncomp

        D_lim1 = -10.
        D_lim2 = -15.
        newtonRaphson = 1
        ncomp= frnw_g(j,1)
        S_ncomp = (-z(ncomp,j)+z(ncomp-1,j))/dx(ncomp-1,j)
        elevTable = xsec_tab(1,:,ncomp,j)
        areaTable = xsec_tab(2,:,ncomp,j)
        topwTable = xsec_tab(6,:,ncomp,j)
        depthCalOk(ncomp) = 1

        call r_interpol(elevTable,areaTable,nel,newY(ncomp,j),newArea(ncomp,j))
        if (newArea(ncomp,j) .eq. -9999) then
!            print*, 'At j = ',j,', i = ',ncomp, 'time =',t, 'interpolation of newArea was not possible'
!            stop
        end if

        call r_interpol(elevTable,topwTable,nel,newY(ncomp,j),bo(ncomp,j))
        if (bo(ncomp,j) .eq. -9999) then
!            print*, 'At j = ',j,', i = ',ncomp, 'time =',t, 'interpolation of bo was not possible'
!            stop
        end if

        do i=ncomp,1,-1
            currentQ = qp(i,j)
            q_sk_multi=1.0
            !* Calculating : read all attributes from tab file
            elevTable = xsec_tab(1,:,i,j)
            convTable = xsec_tab(5,:,i,j)
            areaTable = xsec_tab(2,:,i,j)
            pereTable = xsec_tab(3,:,i,j)
            topwTable = xsec_tab(6,:,i,j)
            skkkTable = xsec_tab(11,:,i,j)
            !* interpolate the cross section attributes based on water elevation
            xt=newY(i,j)
            currentSquareDepth=(elevTable-z(i,j))**2.

            call r_interpol(currentSquareDepth,convTable,nel,(newY(i,j)-z(i,j))**2.0,co(i))

            if (co(i) .eq. -9999) then
!                print*, 'At j = ',j,', i = ',i, 'time =',t, 'interpolation of conveyence was not possible, wl', &
!                newY(i,j), 'z',z(i,j),'previous wl',newY(i+1,j), 'previous z',z(i+1,j), 'dimensionless_D(i,j)', &
!                dimensionless_D(i,j)
!                stop
            end if
            co(i) =q_sk_multi * co(i)

            call r_interpol(elevTable,areaTable,nel,xt,newArea(i,j))

            if (newArea(i,j) .eq. -9999) then
!                print*, 'At j = ',j,', i = ',i, 'time =',t, 'interpolation of newArea was not possible'
!                stop
            end if
            call r_interpol(elevTable,pereTable,nel,xt,pere(i,j))
            call r_interpol(elevTable,topwTable,nel,xt,bo(i,j))
            call r_interpol(elevTable,skkkTable,nel,xt,sk(i,j))

            sfi = qp(i,j) * abs(qp(i,j)) / ( co(i)** 2.0 )
            chnWidth = rightBank(i,j)-leftBank(i,j)
            chnWidth = min(chnWidth,bo(i,j))

            if (depthCalOk(i) .eq. 1) then
                celerity2(i)=5.0 / 3.0 * abs(sfi) ** 0.3 * abs(qp(i,j)) ** 0.4 / bo(i,j) ** 0.4 / (1./(sk(i,j)*q_sk_multi)) ** 0.6
                diffusivity2(i) = abs(qp(i,j)) / 2.0 / bo(i,j) / abs(sfi)
                vel = qp(i,j)/newArea(i,j)
                if (celerity2(i) .gt. 3.0*vel) celerity2(i) = vel*3.0
            else
                if (qp(i,j) .lt. 1) then
                    celerity2(i)=0.5
                else
                    celerity2(i)=1.0
                end if
                diffusivity2(i)=diffusivity(i,j)
            end if
            newMassBalance =0
            if (newMassBalance .eq. 1) then
                if (i .gt. 1) then
                    newArea(i-1,j) = oldArea(i-1,j) + oldArea(i,j) - newArea(i,j) - 2.*dtini/dx(i-1,j)*(qp(i,j)-qp(i-1,j))  ! change 20210407
                    elevTable = xsec_tab(1,:,i-1,j)
                    areaTable = xsec_tab(2,:,i-1,j)
                    if ( newArea(i-1,j) .le. 0) then
                        slope = (z(i-1,j)-z(i,j))/dx(i-1,j)
                        call normal_y(i-1, j, q_sk_multi, slope, qp(i-1,j), newY(i-1,j), temp, newArea(i-1,j), temp)
                        currentRoutingNormal(i-1,j) = 1
                    else
                        call r_interpol(areaTable,elevTable,nel,newArea(i-1,j),newY(i-1,j))
                        currentRoutingNormal(i-1,j) = 0
                    end if
                end if
            else
                if (i .gt. 1) then
                    !! If routing method is changed just a few time steps ago, we maintain the same routing to avoid oscillation
                    if ( (routingNotChanged(i-1,j) .lt. minNotSwitchRouting2) .and. (currentRoutingNormal(i-1,j) .lt. 3) ) then
                        if (currentRoutingNormal(i-1,j) .eq. 0) then
                            newY(i-1,j) = newY(i,j) + sfi * dx(i-1,j)
                        else if (currentRoutingNormal(i-1,j) .eq. 1) then
                            slope = (z(i-1,j)-z(i,j))/dx(i-1,j)
                            if (slope .le. 0.0001) slope = 0.0001
                            q_sk_multi= 1.0
                            ! applying normal depth to all the nodes
                            call normal_crit_y(i-1, j, q_sk_multi, slope, qp(i-1,j), newY(i-1,j), temp, newArea(i-1,j), temp)
                        end if
                    else
                        !! If DSP: D is below 1.0, we switch to partial diffusive routing
                        if (dimensionless_D(i-1,j) .lt. D_lim1) then
                            slope = (z(i-1,j)-z(i,j))/dx(i-1,j)
                            if (slope .le. 0.0001) slope = 0.0001
                                q_sk_multi= 1.0
                                ! applying normal depth to all the nodes
                                call normal_crit_y(i-1, j, q_sk_multi, slope, qp(i-1,j), newY(i-1,j), temp, newArea(i-1,j), temp)
                                 ! Book-keeping: changing from full diffusive to partial diffusive
                                if ( currentRoutingNormal(i-1,j) .ne. 1 ) routingNotChanged(i-1,j) = 0
                                currentRoutingNormal(i-1,j) = 1
                            !! If DSP: D is not below 1.0, we switch to full diffusive routing
                            elseif ( (dimensionless_D(i-1,j) .ge. D_lim1) .and. (dimensionless_D(i-1,j) .lt. D_lim2) ) then
!                                print*, 'partial diffusive at j', j, 'i-1',i-1
                                slope = (z(i-1,j)-z(i,j))/dx(i-1,j)
                                if (slope .le. 0.0001) slope = 0.0001
                                q_sk_multi=1.0
                                ! applying normal depth to all the nodes
                                call normal_crit_y(i-1, j, q_sk_multi, slope, qp(i-1,j), stg1, temp, newArea(i-1,j), temp)

                                stg2 = newY(i,j) + sfi * dx(i-1,j)
                                newY(i-1,j) = ( stg2 * (dimensionless_D(i-1,j) - D_lim1) + &
                                                stg1 * (D_lim2 - dimensionless_D(i-1,j)) ) / (D_lim2 - D_lim1)
                                 ! Book-keeping: changing from full diffusive to partial diffusive
                                if ( currentRoutingNormal(i-1,j) .ne. 3 ) routingNotChanged(i-1,j) = 0
                                currentRoutingNormal(i-1,j) = 3
                            else
                                slope = (z(i-1,j)-z(i,j))/dx(i-1,j)
                                depthYi = newY(i,j) - z(i,j)
                                tempDepthi_1 = oldY(i-1,j)-z(i-1,j)
                                elevTable_1 = xsec_tab(1,:,i-1,j)
                                areaTable_1 = xsec_tab(2,:,i-1,j)
                                pereTable_1 = xsec_tab(3,:,i-1,j)
                                rediTable_1 = xsec_tab(4,:,i-1,j)
                                convTable_1 = xsec_tab(5,:,i-1,j)
                                topwTable_1 = xsec_tab(6,:,i-1,j)
                                dKdATable_1 = xsec_tab(9,:,i-1,j)       ! change Nazmul 20210601
                                skkkTable_1 = xsec_tab(11,:,i-1,j)
                                currentSquaredDepth_1=(elevTable_1-z(i-1,j))**2.0
                                toll = 1.0
                                iii = 0
                                ! Applying Newton\96Raphson method
                                if (newtonRaphson .eq. 1) then
                                    do while ( abs(toll) .gt. 0.001)
                                        iii = iii +1
                                        tempY_1 = tempDepthi_1 + z(i-1,j)

                                        call r_interpol(currentSquaredDepth_1,convTable_1,nel,(tempDepthi_1)**2.0,tempCo_1)

                                        temp_q_sk_multi_1=1.0
                                        tempCo_1 = tempCo_1 * temp_q_sk_multi_1

                                        call r_interpol(elevTable_1,areaTable_1,nel,tempY_1,tempArea_1)
                                        call r_interpol(elevTable_1,pereTable_1,nel,tempY_1,tempPere_1)
                                        call r_interpol(elevTable_1,rediTable_1,nel,tempY_1,tempRadi_1)
                                        call r_interpol(elevTable_1,topwTable_1,nel,tempY_1,tempbo_1)
                                        call r_interpol(elevTable_1,dKdATable_1,nel,tempY_1,tempdKdA_1)     ! change Nazmul 20210601

                                        tempsfi_1 = qp(i-1,j) * abs(qp(i-1,j)) / ( tempCo_1** 2.0 )!

                                        ffy = tempDepthi_1 - depthYi + dx(i-1,j)*slope - 0.5*dx(i-1,j)*(sfi + tempsfi_1)

                                        dkda = tempdKdA_1                                                   ! change Nazmul 20210601

                                        ffprime = 1 + dx(i-1,j) * tempbo_1 *  qp(i-1,j) * abs(qp(i-1,j)) / (tempCo_1 ** 3.0) * dkda

                                        tempDepthi_1_new = tempDepthi_1 - ffy / ffprime

                                        tempDepthi_1_new = max(tempDepthi_1_new,0.005)

                                        toll = abs(tempDepthi_1_new - tempDepthi_1)

                                        ! Change Nazmul 20210601
                                        if(iii .gt. 30)then
!                                            print*, 'Warning: Depth iteration reached maximum trial at j=', j, 'i=', i-1 , &
!                                            'and',i,'depths are',tempDepthi_1, tempDepthi_1_new, 'slope=', slope, &
!                                            'dx=', dx(i-1,j), 'depth at d/s',depthYi, 'Q-s are', qp(i-1,j), qp(i,j)
                                            depthCalOk(i-1) = 0
                                            EXIT
                                        endif

                                        tempDepthi_1 = tempDepthi_1_new
                                        depthCalOk(i-1) = 1
                                    end do
                                end if

                                ! Applying mid point bisection
                                if (newtonRaphson .eq. 0) then
                                    tempY_1 = elevTable_1(2)
                                    tempY_2 = depthYi * 3. + z(i-1,j)
                                    tempY_3 = (tempY_1 + tempY_2) / 2.
                                    do while ( abs(toll) .gt. 0.001)
                                        iii = iii +1

                                        call r_interpol(currentSquaredDepth_1,convTable_1,nel,(tempY_1-z(i-1,j))**2.0,tempCo_1)
                                        call r_interpol(currentSquaredDepth_1,convTable_1,nel,(tempY_2-z(i-1,j))**2.0,tempCo_2)
                                        call r_interpol(currentSquaredDepth_1,convTable_1,nel,(tempY_3-z(i-1,j))**2.0,tempCo_3)

                                        temp_q_sk_multi_1=1.0
                                        tempCo_1 = tempCo_1 * temp_q_sk_multi_1
                                        tempCo_2 = tempCo_2 * temp_q_sk_multi_1
                                        tempCo_3 = tempCo_3 * temp_q_sk_multi_1

                                        tempsfi_1 = qp(i-1,j) * abs(qp(i-1,j)) / ( tempCo_1** 2.0 )
                                        tempsfi_2 = qp(i-1,j) * abs(qp(i-1,j)) / ( tempCo_2** 2.0 )
                                        tempsfi_3 = qp(i-1,j) * abs(qp(i-1,j)) / ( tempCo_3** 2.0 )

                                        ffy_1 = (tempY_1-z(i-1,j)) - depthYi + dx(i-1,j)*slope - 0.5*dx(i-1,j)*(sfi + tempsfi_1)
                                        ffy_2 = (tempY_2-z(i-1,j)) - depthYi + dx(i-1,j)*slope - 0.5*dx(i-1,j)*(sfi + tempsfi_2)
                                        ffy_3 = (tempY_3-z(i-1,j)) - depthYi + dx(i-1,j)*slope - 0.5*dx(i-1,j)*(sfi + tempsfi_3)

                                        if ((ffy_1 * ffy_2) .gt. 0.) then
                                            tempY_2 = (tempY_2 - z(i-1,j)) * 2.0 + z(i-1,j)
                                        elseif ((ffy_1 * ffy_3) .le. 0.) then
                                            tempY_2 = tempY_3
                                        elseif ((ffy_2 * ffy_3) .le. 0.) then
                                            tempY_1 = tempY_3
                                        end if
                                        tempY_3 = (tempY_1 + tempY_2) / 2.0
                                        toll = tempY_2 - tempY_1
                                        tempDepthi_1 = tempY_3 - z(i-1,j)
                                        depthCalOk(i-1) = 1
                                    end do
                                end if
                                newY(i-1,j) = tempDepthi_1 + z(i-1,j)

                                if (newY(i-1,j) .gt. 10.0**5.) newY(i-1,j) = 10.0**5.
                                 ! Book-keeping: changing from partial diffusive to full diffusive
                                if ( currentRoutingNormal(i-1,j) .ne. 0 ) routingNotChanged(i-1,j) = 0
                                currentRoutingNormal(i-1,j) = 0
                            end if
                        end if

                        if (newY(i-1,j)-z(i-1,j) .le. 0.) then
                            !print*, 'depth is negative at time=,', t,'j= ', j,'i=',i-1,'newY=',(newY(jj,j),jj=1,ncomp)
!                            print*, 'dimensionless_D',(dimensionless_D(jj,j),jj=1,ncomp)
!                            print*, 'newQ',(newQ(jj,j),jj=1,ncomp)
!                            print*, 'Bed',(z(jj,j),jj=1,ncomp)
!                            print*, 'dx',(dx(jj,j),jj=1,ncomp-1)
                            !pause 777
                        end if
                    end if
                end if

                ! Book-keeping: Counting the number as for how many time steps the routing method is unchanged
                if (i.gt.1) then
                    routingNotChanged(i-1,j) = routingNotChanged(i-1,j) + 1
                endif
            end do

            celerity(1:ncomp,j) =  sum(celerity2(1:ncomp)) / ncomp  ! change Nazmul 20210601
            if (celerity(1,j) .lt. 0.5) celerity(1:ncomp,j) = 0.5
            diffusivity(1:ncomp,j)=sum(diffusivity2(1:ncomp)) / ncomp
            do i = 1, ncomp
                if (diffusivity(i,j) .gt. 1000.) diffusivity(i,j) = 1000. !!! Test
                if (diffusivity(i,j) .lt. 50.) diffusivity(i,j) = 50. !!! Test
            end do
    end subroutine mesh_diffusive_backward
    !**-----------------------------------------------------------------------------------------
    !*      Create lookup tables at each node storing computed values of channel geometries
    !*      such as area and conveyance and normal/critical depth for possible ranges of
    !*      water depth.
    !
    !**-----------------------------------------------------------------------------------------
    subroutine readXsection(k,lftBnkMann,rmanning_main,rgtBnkMann,leftBnkX_given,rghtBnkX_given,timesDepth,num_reach,&
                            z_ar_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g )
        implicit none
        save

        integer, intent(in) :: k, num_reach
        doubleprecision, intent(in) :: rmanning_main,lftBnkMann,rgtBnkMann,leftBnkX_given,rghtBnkX_given, timesDepth
        double precision, dimension(mxncomp, nlinks), intent(in) :: z_ar_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g
        doubleprecision, dimension(:), allocatable :: xcs, ycs
        doubleprecision, dimension(:,:), allocatable :: el1, a1, peri1, redi1
        doubleprecision, dimension(:), allocatable :: redi1All
        doubleprecision, dimension(:,:), allocatable :: conv1, tpW1, diffArea, newI1, diffPere
        doubleprecision, dimension(:), allocatable :: newdPdA, diffAreaAll, diffPereAll, newdKdA       ! change Nazmul 20210601
        doubleprecision, dimension(:), allocatable :: compoundSKK, elev
        integer, dimension(:), allocatable :: i_start, i_end, totalNodes
        doubleprecision, dimension(:,:), allocatable :: allXcs, allYcs
        integer :: i_area, i_find, i, j, jj, num
        doubleprecision :: el_min, el_max, el_range, el_incr, el_now, x1, y1, x2, y2, x_start, x_end, waterElev, leftBnkX,rghtBnkX
        doubleprecision :: f2m, cal_area, cal_peri, cal_topW,  diffAreaCenter
        doubleprecision :: compoundMann, el_min_1
        integer:: i1, i2
        doubleprecision :: leftBnkY, rghtBnkY,rmanning
        integer:: mainChanStrt, mainChanEnd,  kkk, startFound, endFound
        doubleprecision :: hbf

        allocate (el1(nel,3),a1(nel,3),peri1(nel,3),redi1(nel,3),redi1All(nel))
        allocate (conv1(nel,3), tpW1(nel,3), diffArea(nel,3), newI1(nel,3), diffPere(nel,3))
        allocate (newdPdA(nel), diffAreaAll(nel), diffPereAll(nel), newdKdA(nel))       ! change Nazmul 20210601
        allocate (compoundSKK(nel), elev(nel))
        allocate (i_start(nel), i_end(nel))
        allocate (totalNodes(3))


        leftBnkX=leftBnkX_given
        rghtBnkX=rghtBnkX_given
        startFound = 0
        endFound = 0
        !* channel geometry at a given segment
        z_g = z_ar_g(k, num_reach)
        bo_g= bo_ar_g(k, num_reach)
        traps_g=traps_ar_g(k, num_reach)
        tw_g= tw_ar_g(k, num_reach)
        twcc_g= twcc_ar_g(k, num_reach)
        hbf= (tw_g-bo_g)/(2.0*traps_g) !* bankfull depth
        maxTableLength=8
        f2m=1.0
        allocate (xcs(maxTableLength), ycs(maxTableLength))
        allocate (allXcs(maxTableLength,3), allYcs(maxTableLength,3))
        do i=1, maxTableLength
            !* channel x-section vertices at a given segment
            if (i==1) then
                x1=0.0; y1=z_g + timesDepth*hbf
            elseif (i==2) then
                x1=0.0; y1=z_g + hbf
            elseif (i==3) then
                x1=(twcc_g-tw_g)/2.0; y1= z_g + hbf
            elseif (i==4) then
                x1= xcs(3) + traps_g*hbf; y1= z_g
            elseif (i==5) then
                x1= xcs(4) + bo_g; y1= z_g
            elseif (i==6) then
                x1= xcs(5) + traps_g*hbf; y1= z_g + hbf
            elseif (i==7) then
                x1= twcc_g; y1= z_g + hbf
            elseif (i==8) then
                x1= xcs(7); y1= z_g + timesDepth*hbf
            endif

            xcs(i)=x1*f2m
            ycs(i)=y1*f2m
            if ((xcs(i) .ge. leftBnkX) .and. (startFound .eq. 0)) then
                mainChanStrt = i-1
                startFound = 1
            end if
            if ((xcs(i) .ge. rghtBnkX) .and. (endFound .eq. 0)) then
                mainChanEnd = i-1
                endFound = 1
            end if
        enddo
        mainChanStrt=3
        mainChanEnd=6
        num=i

        if (leftBnkX .lt. minval(xcs(2:num-1))) leftBnkX = minval(xcs(2:num-1))
        if (rghtBnkX .gt. maxval(xcs(2:num-1))) rghtBnkX = maxval(xcs(2:num-1))

        leftBnkY = ycs(mainChanStrt)+(leftBnkX-xcs(mainChanStrt))/&
          (xcs(mainChanStrt+1)-xcs(mainChanStrt))*(ycs(mainChanStrt+1)-ycs(mainChanStrt))
        rghtBnkY = ycs(mainChanEnd)+(rghtBnkX-xcs(mainChanEnd))/&
          (xcs(mainChanEnd+1)-xcs(mainChanEnd))*(ycs(mainChanEnd+1)-ycs(mainChanEnd))
        el_min=99999.
        el_max=-99999.
        do i=2,num-1
            if(ycs(i).lt.el_min)el_min=ycs(i)
            if(ycs(i).gt.el_max)el_max=ycs(i)
        enddo
        el_range=(el_max-el_min)*2.0 ! change Nazmul 20210601

        do i=1, 3
            allXcs(i+1,1)=xcs(i) !x1*f2m
            allYcs(i+1,1)=ycs(i) !y1*f2m
        enddo
        allXcs(1,1)=xcs(1)
        allYcs(1,1)=el_min+el_range+1.
        allXcs(mainChanStrt+2,1)=xcs(3)
        allYcs(mainChanStrt+2,1)=el_min+el_range+1.

        do i=3,4
            allXcs(i-1,2)=xcs(i) !x1*f2m
            allYcs(i-1,2)=ycs(i) !y1*f2m
        enddo

        do i=5,6
            allXcs(i,2)=xcs(i) !x1*f2m
            allYcs(i,2)=ycs(i) !y1*f2m
        enddo
        allXcs(1,2)=xcs(3)
        allYcs(1,2)=el_min+el_range+1.
        allXcs(7,2)=xcs(6)
        allYcs(7,2)=el_min+el_range+1.

        do i=6,8
            allXcs(i-4,3)=xcs(i) !x1*f2m
            allYcs(i-4,3)=ycs(i) !y1*f2m
        enddo
        allXcs(1,3) = allXcs(2,3)
        allYcs(1,3) = el_min+el_range+1.
        i=5
        allXcs(i,3) = allXcs(i-1,3)
        allYcs(i,3) = el_min+el_range+1.

        totalNodes(1) = 5
        totalNodes(2) = 7
        totalNodes(3) = 5

        allXcs(4,2) = (allXcs(3,2)+allXcs(5,2))/2.0
        allYcs(4,2) = allYcs(3,2) - 0.01

        el_min_1 = el_min
        el_min = allYcs(4,2)    ! change Nazmul 20210601 ! corrected

        elev(1) = el_min
        elev(2) = el_min + 0.01/4.
        elev(3) = el_min + 0.01/4.*2.
        elev(4) = el_min + 0.01/4.*3.
        elev(5) = el_min + 0.01

        el_incr=el_range/real(nel-6.0)

        do kkk = 6,nel
            elev(kkk) = elev(5)+el_incr * (kkk-5)
        end do

        xcs = 0.
        ycs = 0.
        newI1=0.0 !Hu changed
        do kkk=1,3
            num = totalNodes(kkk)
            xcs(1:num) = allXcs(1:num,kkk)
            ycs(1:num) = allYcs(1:num,kkk)
            if (kkk .eq. 1) rmanning = lftBnkMann
            if (kkk .eq. 2) rmanning = rmanning_main
            if (kkk .eq. 3) rmanning = rgtBnkMann
            do j=1,nel
                el_now = elev(j)
                if(abs(el_now - el_min) < TOLERANCE) then
                    el_now=el_now+0.00001
                end if
                i_start(1)=-999
                i_end(1)=-999
                i_area=0
                i_find=0
                do i=1,num-1
                    y1=ycs(i)
                    y2=ycs(i+1)
                    if(el_now.le.y1 .and. el_now.gt.y2 .and. i_find.eq.0)then
                        i_find=1
                        i_area=i_area+1
                        i_start(i_area)=i
                    endif
                    if(el_now.gt.y1 .and. el_now.le.y2 .and. i_find.eq.1)then
                        i_find=0
                        i_end(i_area)=i
                    endif
                enddo

                cal_area=0.
                cal_peri=0.
                cal_topW=0.

                do i=1,i_area
                    x1=xcs(i_start(i))
                    x2=xcs(i_start(i)+1)
                    y1=ycs(i_start(i))
                    y2=ycs(i_start(i)+1)
                    if(y1.eq.y2)then
                        x_start=x1
                    else
                        x_start=x1+(el_now-y1)/(y2-y1)*(x2-x1)
                    endif
                    write(11,*)x_start,el_now

                    x1=xcs(i_end(i))
                    x2=xcs(i_end(i)+1)
                    y1=ycs(i_end(i))
                    y2=ycs(i_end(i)+1)

                    if(y1.eq.y2)then
                      x_end=x1
                    else
                      x_end=x1+(el_now-y1)/(y2-y1)*(x2-x1)
                    endif

                    cal_topW=x_end-x_start+cal_topW

                    write(11,*)x_end,el_now
                    write(11,*)'NaN NaN'

                    i1=i_start(i)
                    i2=i_end(i)

                    cal_area = cal_area    &
                             +cal_tri_area(el_now,x_start,xcs(i1+1),ycs(i1+1))    &
                             +cal_multi_area(el_now,xcs,ycs,maxTableLength,i1+1,i2)    &
                             +cal_tri_area(el_now,x_end,xcs(i2),ycs(i2))
                    cal_peri = cal_peri    &
                            +cal_dist(x_start,el_now,xcs(i1+1),ycs(i1+1))    &
                            +cal_perimeter(xcs,ycs,maxTableLength,i1+1,i2)    &
                            +cal_dist(x_end,el_now,xcs(i2),ycs(i2))
                    if(i1.eq.1)cal_peri=cal_peri    &
                             -cal_dist(x_start,el_now,xcs(i1+1),ycs(i1+1))
                    if(i2.eq.(num-1))cal_peri=cal_peri    &
                             -cal_dist(x_end,el_now,xcs(i2),ycs(i2))

                enddo

                el1(j,kkk)=el_now
                a1(j,kkk)=cal_area
                peri1(j,kkk)=cal_peri
                redi1(j,kkk)=a1(j,kkk)/peri1(j,kkk)

                conv1(j,kkk)=1./rmanning*a1(j,kkk)*(redi1(j,kkk))**(2./3.)
                if (peri1(j,kkk) .le. TOLERANCE) then
                    redi1(j,kkk) =0.0; conv1(j,kkk)=0.0
                endif
                tpW1(j,kkk)=cal_topW

                if(j.eq.1) then !* Dongha added
                    diffArea(j,kkk)=a1(j,kkk) !* Dongha added
                    diffPere(j,kkk)=peri1(j,kkk) !* Dongha added
                else
                    if (el_now .le. minval(ycs(1:num))) then
                      diffArea(j,kkk)=a1(j,kkk)
                      diffPere(j,kkk)=peri1(j,kkk)
                    else
                      diffArea(j,kkk)=a1(j,kkk)-a1(j-1,kkk)
                      diffPere(j,kkk)=peri1(j,kkk)-peri1(j-1,kkk)
                    endif
                endif

                waterElev=el1(j,kkk)
                do jj=2,j
                  diffAreaCenter=el1(jj,kkk)-(el1(jj,kkk)-el1(jj-1,kkk))*0.5
                  newI1(j,kkk)=newI1(j,kkk)+diffArea(jj,kkk)*(waterElev-diffAreaCenter)
                enddo
            end do
        end do

        do j = 1,nel
            el_now=el1(j,1)
            if (j .eq. 1) then
                newdPdA(j) = sum(peri1(j,:)) / sum(a1(j,:))
                newdKdA(j) = sum(conv1(j,:)) / sum(a1(j,:))     ! change Nazmul 20210601
            else
                newdPdA(j)= (sum(peri1(j,:)) - sum(peri1(j-1,:))) / (sum(a1(j,:)) - sum(a1(j-1,:)))
                newdKdA(j)= (sum(conv1(j,:)) - sum(conv1(j-1,:))) / (sum(a1(j,:)) - sum(a1(j-1,:)))
            end if

            compoundMann = sqrt((abs(peri1(j,1))*lftBnkMann ** 2. + abs(peri1(j,2))*rmanning_main ** 2.+&
             abs(peri1(j,3))*rgtBnkMann ** 2.) / (abs(peri1(j,1))+abs(peri1(j,2))+abs(peri1(j,3))))
            compoundSKK(j) = 1. / compoundMann

            redi1All(j)=sum(a1(j,:)) /sum(peri1(j,:))
            xsec_tab(1,j,k,num_reach) = el1(j,1)
            xsec_tab(2,j,k,num_reach) = sum(a1(j,:))
            xsec_tab(3,j,k,num_reach) = sum(peri1(j,:))
            xsec_tab(4,j,k,num_reach) = redi1All(j)
            xsec_tab(5,j,k,num_reach) = sum(conv1(j,:))
            xsec_tab(6,j,k,num_reach) = abs(tpW1(j,1))+abs(tpW1(j,2))+abs(tpW1(j,3))
            xsec_tab(7,j,k,num_reach) = sum(newI1(j,:))
            xsec_tab(8,j,k,num_reach) = newdPdA(j)
            xsec_tab(9,j,k,num_reach) = newdKdA(j)
            xsec_tab(11,j,k,num_reach) = compoundSKK(j)
        end do
        z(k,num_reach)= el_min

        deallocate (el1, a1, peri1, redi1, redi1All)
        deallocate (conv1, tpW1, diffArea, newI1, diffPere)
        deallocate (newdPdA, diffAreaAll, diffPereAll, newdKdA)       ! change Nazmul 20210601
        deallocate (compoundSKK, elev)
        deallocate (i_start, i_end)
        deallocate (totalNodes)
        deallocate (xcs, ycs)
        deallocate (allXcs, allYcs)

        contains
            !**----------------------------------------
            !*      calculate area of triangle
            !**----------------------------------------
            double precision function cal_tri_area(el,x0,x1,y1)
                implicit none
                  doubleprecision, intent(in) :: el,x0,x1,y1
                  cal_tri_area=abs(0.5*(x1-x0)*(el-y1))
                  return
            end function cal_tri_area
            !**----------------------------------------
            !*      calculate area of trapezoid
            !**----------------------------------------
            double precision function cal_trap_area(el,x1,y1,x2,y2)
                implicit none
                doubleprecision, intent(in) :: el,x1,y1,x2,y2
                cal_trap_area=abs(0.5*(x2-x1)*(el-y1+el-y2))
                return
            end function cal_trap_area
            !**--------------------------------------------
            !*     calculate sum of areas of trapezoids
            !**--------------------------------------------
            double precision function cal_multi_area(el,xx,yy,n,i1,i2)
                implicit none
                integer, intent(in) :: n,i1,i2
                doubleprecision, intent(in) :: el
                doubleprecision, intent(in) :: xx(n),yy(n)
                integer :: i
                doubleprecision :: area, x1, x2, y1, y2
                area=0
                do i=i1,i2-1
                    x1=xx(i)
                    y1=yy(i)
                    x2=xx(i+1)
                    y2=yy(i+1)
                    area=area+cal_trap_area(el,x1,y1,x2,y2)
                enddo
                cal_multi_area=area
                return
            endfunction cal_multi_area
            !**--------------------------------------------
            !*     calculate distance of two vertices
            !**--------------------------------------------
            double precision function cal_dist(x1,y1,x2,y2)
                implicit none
                doubleprecision, intent(in) :: x1,y1,x2,y2
                cal_dist=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+1.e-32)
                return
            end function cal_dist
            !**--------------------------------------------
            !*     calculate wetted perimeter
            !**--------------------------------------------
            double precision function cal_perimeter(xx,yy,n,i1,i2)
                implicit none
                integer, intent(in) :: n,i1,i2
                doubleprecision, intent(in) :: xx(n),yy(n)
                integer :: i
                doubleprecision :: p, x1, x2, y1, y2
                p=0.
                do i=i1,i2-1
                    x1=xx(i)
                    y1=yy(i)
                    x2=xx(i+1)
                    y2=yy(i+1)
                    p=p+cal_dist(x1,y1,x2,y2)
                enddo
                cal_perimeter=p
                return
            endfunction cal_perimeter
    endsubroutine readXsection
    !*---------------------------------------------------
    !*      interpolation with given arrays
    !
    !*---------------------------------------------------
    subroutine r_interpol(x,y,kk,xrt,yt)
        implicit none
        integer, intent(in) :: kk
        doubleprecision, intent(in) :: xrt, x(kk), y(kk)
        doubleprecision, intent(out) :: yt
        integer :: k

        if (xrt.le.maxval(x) .and. xrt.ge.minval(x)) then
            do k=1,kk-1
                if((x(k)-xrt)*(x(k+1)-xrt).le.0)then

                    yt=(xrt-x(k))/(x(k+1)-x(k))*(y(k+1)-y(k))+y(k)

                    EXIT
                endif
            end do
        else if (xrt.ge.maxval(x)) then
            !print*, xrt, ' is above the user defined limit'
            yt=(xrt-x(kk-1))/(x(kk)-x(kk-1))*(y(kk)-y(kk-1))+y(kk-1) ! extrapolation

        else
            print*, xrt, ' is below the user defined limit'
            yt = -9999.0

            print*, 'maxval(x)= ', maxval(x), 'and minval(x)=', minval(x),'so',  xrt, ' is not within the limit'
            print*, 'kk', kk
            print*, 'x', (x(k), k=1, kk)
            print*, 'y', (y(k), k=1, kk)
        end if
    end subroutine r_interpol
    !*----------------------------------------------------------------
    !*     compute normal/critical depth/area using lookup tables
    !
    !*----------------------------------------------------------------
    subroutine normal_crit_y(i, j, q_sk_multi, So, dsc, y_norm, y_crit, area_n, area_c)
        implicit none
        integer, intent(in) :: i, j
        doubleprecision, intent(in) :: q_sk_multi, So, dsc
        doubleprecision, intent(out) :: y_norm, y_crit, area_n, area_c
        doubleprecision :: area_0, width_0, errorY, pere_0, hydR_0, skk_0
        integer :: trapnm_app, recnm_app, iter

        elevTable = xsec_tab(1,:,i,j)
        areaTable = xsec_tab(2,:,i,j)
        pereTable = xsec_tab(3,:,i,j)
        convTable = xsec_tab(5,:,i,j)
        topwTable = xsec_tab(6,:,i,j)
        call r_interpol(convTable,areaTable,nel,dsc/sqrt(So),area_n)
        call r_interpol(convTable,elevTable,nel,dsc/sqrt(So),y_norm)
        call r_interpol(elevTable,areaTable,nel,oldY(i,j),area_0) ! initial estimate
        call r_interpol(elevTable,topwTable,nel,oldY(i,j),width_0) ! initial estimate

        area_c=area_0
        errorY = 100.
        !pause
        do while (errorY .gt. 0.0001)
            area_c = (dsc * dsc * width_0 / grav) ** (1./3.)
            errorY = abs(area_c - area_0)
            call r_interpol(areaTable,topwTable,nel,area_c, width_0)
            area_0 = area_c
        enddo

        call r_interpol(areaTable,elevTable,nel,area_c,y_crit)
        if (y_norm .eq. -9999) then
!            print*, 'At j = ',j,', i = ',i, 'interpolation of y_norm in calculating normal area was not possible, Q', &
!            dsc,'slope',So !,'lateralFlow', lateralFlow(1:nx1(j),j)
!            stop
        end if
    end subroutine normal_crit_y
    !*--------------------------------------------------
    !*     compute normal depth using lookup tables
    !
    !*--------------------------------------------------
    subroutine normal_y(i, j, q_sk_multi, So, dsc, y_norm, y_crit, area_n, area_c)
        implicit none
        integer, intent(in) :: i, j
        doubleprecision, intent(in) :: q_sk_multi, So, dsc
        doubleprecision, intent(out) :: y_norm, y_crit, area_n, area_c
        doubleprecision :: area_0, width_0, errorY, hydR_0,skk_0!, fro
        integer :: trapnm_app, recnm_app, iter

        elevTable = xsec_tab(1,:,i,j)
        areaTable = xsec_tab(2,:,i,j)
        rediTable = xsec_tab(4,:,i,j)
        topwTable = xsec_tab(6,:,i,j)
        skkkTable = xsec_tab(11,:,i,j)

        call r_interpol(elevTable,areaTable,nel,oldY(i,j),area_0) ! initial estimate

        errorY = 100.
        do while (errorY .gt. 0.00001)
            call r_interpol(areaTable,rediTable,nel,area_0,hydR_0)
            call r_interpol(areaTable,skkkTable,nel,area_0,skk_0)
            area_n = dsc/skk_0/q_sk_multi/ hydR_0 ** (2./3.) / sqrt(So)
            errorY = abs(area_n - area_0) / area_n
            area_0 = area_n
        enddo
        call r_interpol(areaTable,elevTable,nel,area_0,y_norm)
        y_crit = -9999.
        area_c = -9999.
    end subroutine normal_y
    !*--------------------------------------------------
    !*                 Linear Interpolation
    !
    !*--------------------------------------------------
    double precision function LInterpol(x1,y1,x2,y2,x)
        implicit none
        doubleprecision, intent(in) :: x1, y1, x2, y2, x
        !* interpolate y for the given x
        LInterpol= (y2-y1)/(x2-x1)*(x-x1)+y1
    end function LInterpol
    !*--------------------------------------------
    !           Interpolate any value
    !
    !*--------------------------------------------
    double precision function intp_y(nrow, xarr, yarr, x)
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
        doubleprecision, dimension(:), intent(in) :: xx
        doubleprecision, intent(in) :: x
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
endmodule mdiffv2





























