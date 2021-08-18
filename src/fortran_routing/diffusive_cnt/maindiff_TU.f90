module diffusive
    !*-------------------------------------------------------------------------------------------------
    !*       A diffusive model developed by Tulane University (Prof.Ehab Meselhe and
    !*       Dr.Md Nazmul Azim Beg) and integrated into a Python routing framework by Dong Ha Kim at
    !*       the National Water Center. Basically, the partial differential equation of diffusive
    !*       wave is numerically solved using Crank-Nicoloson and Hermite Interpolation techniques.
    !*       Water depth computation produces either normal depth or diffusive depth, the selection
    !*       of which is determined mainly by dimensionless diffusion coefficient.
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
    double precision :: dtini, dxini, cfl, minDx, maxCelerity,  theta, min_Q					! change Nazmul CNT
    double precision :: frus2, minNotSwitchRouting, minNotSwitchRouting2

    double precision, dimension(:), allocatable :: area, depth, co, froud, courant
    double precision, dimension(:,:), allocatable :: bo, dx
	!**arrays for branching channel application
    double precision, dimension(:,:), allocatable :: areap, z, sk								! change Nazmul CNT
	double precision, dimension(:,:,:), allocatable :: qp										! change Nazmul CNT
    double precision, dimension(:,:), allocatable :: dqp, dap, dqc, dac
    double precision, dimension(:,:), allocatable :: celerity, velocity, diffusivity, qpx		! change Nazmul CNT
    double precision, dimension(:), allocatable :: eei, ffi, exi, fxi, qcx, diffusivity2, celerity2
    ! change for unsteady flow
    double precision, dimension(:,:), allocatable :: pere, oldQ, oldArea, newArea, oldY, newY	! change Nazmul CNT
	double precision, dimension(:,:,:), allocatable :: newQ, lateralFlow, added_Q				! change Nazmul CNT
	double precision, dimension(:,:), allocatable :: ini_q_repeat, ini_E, ini_F					! change Nazmul CNT
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
    double precision :: dmyt, dmyi, dmyj

contains
    !*--------------------------------------------------------------------------------
    !*       Route network by using Python-to-Fortran network traversal map together
    !*          with diffusive routing engines
    !
    !*--------------------------------------------------------------------------------
	
	! TO DO:
	! * dtini_g == saveinterval_g == saveinterval_ev_g. Remove this argument redundancy
	
    subroutine diffnw(dtini_g, t0_g, tfin_g, saveinterval_g, saveinterval_ev_g, dt_ql_g, dt_ub_g, dt_db_g, &
                        nts_ql_g, nts_ub_g, nts_db_g, &
                        mxncomp_g, nrch_g, z_ar_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g, &
                        mann_ar_g, manncc_ar_g, so_ar_g, dx_ar_g, iniq, &
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
        double precision, dimension(mxncomp_g, nrch_g), intent(in) :: mann_ar_g, manncc_ar_g, dx_ar_g, iniq
        double precision, dimension(mxncomp_g, nrch_g, nhincr_m_g), intent(in) :: ufhlt_m_g,  ufqlt_m_g
        double precision, dimension(mxncomp_g, nrch_g, nhincr_f_g), intent(in) :: ufhlt_f_g, ufqlt_f_g

        double precision, dimension(nrch_g, frnw_col), intent(in) :: dfrnw_g !* set frnw_col=10
        double precision, dimension(nts_ql_g, mxncomp_g, nrch_g), intent(in) :: qlat_g
        double precision, dimension(nts_ub_g, nrch_g), intent(in) :: ubcd_g
        double precision, dimension(nts_db_g), intent(in) :: dbcd_g

        double precision, intent(in) :: cfl_g, theta_g, so_llm_g
        integer, intent(in) :: tzeq_flag_g !* 0 for lookup tabale; 1 for using procedures to compute trapz.ch.geo.
        integer, intent(in) :: y_opt_g  !* 1 for normal depth(kinematic); 2 for dept of diffusive wave.

        double precision, dimension(mxncomp_g, nrch_g), intent(inout) :: so_ar_g
		double precision, dimension(ntss_ev_g, mxncomp_g, nrch_g), intent(out) :: q_ev_g, elv_ev_g
        integer :: ncomp

        integer :: i, j, k, ppn, qqn, n, ntim, num_time, igate, pp, boundaryFileMaxEntry, saveFrequency
        integer :: linknb_ds, linknb_us, linknb, nodenb
        double precision :: qnp1_ds, qnp1_us, qsum, y_ds
        double precision :: cour, da, dq, x, saveInterval, width
        double precision :: qn, xt, maxCourant, dtini_given
        double precision :: frds, areasum, yk_ncomp, yav, areak_ncomp, areav, sumOldQ, currentQ, area_ds
        double precision :: arean, areac, hyrdn, hyrdc, perimn, perimc, qcrit, s0ds, timesDepth
        doubleprecision :: latFlowValue, latFlowValue2
        double precision :: t, r_interpol_time, tfin, t1, t2, t0, ini_time !t0 start time
        integer :: tableLength, timestep, kkk, repeatInterval, totalTimeSteps
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
        double precision :: slope, y_norm, area_n, temp

        ! simulation timestep
        dtini=dtini_g
        dtini_given=dtini
		! simulation intial time (hrs)
        t0=t0_g
		! simulation final time (hrs)
        tfin=tfin_g
		! Total number of timesteps in the simulation
		totalTimeSteps = floor((tfin - t0)/dtini*3600)+1 
		! repeat interval (number of timesteps)
		repeatInterval = int(60.*60./dtini_given) ! TO DO, define in terms of a parameter
		! number of timesteps per repeatInterval
        ntim=repeatInterval+5 ! +1 because of boundary effects, at least 1 extra necessary
		num_time=ntim
		! number of reaches in the network
		nlinks=nrch_g
		! network metadata array
        allocate(frnw_g(nlinks,frnw_col))
        frnw_g=dfrnw_g
		! maximum number of segments in a single reach
        mxncomp=mxncomp_g
		!* water depth multiplier used in readXsection
        timesDepth= 4.0 
		! number of rows in channel geometry attribute tables
        nel= 501
		! saving interval (seconds)
        saveInterval=saveInterval_ev_g
		! number of simulation timesteps per save
        saveFrequency=saveInterval / dtini_given
		! maximum number of segments in a single reach - redundant?
        num_points= mxncomp
		! number of reaches in the network - redundant?
        totalChannels= nlinks
		
		
        allocate(area(num_points))
        allocate(bo(num_points,totalChannels))
        allocate(pere(num_points,totalChannels))
        allocate(areap(num_points,totalChannels))
        allocate(qp(num_points,num_time,totalChannels))				! change Nazmul CNT	
        allocate(ini_q_repeat(num_points,totalChannels))			! change Nazmul CNT	
        allocate(ini_E(num_points,totalChannels))					! change Nazmul CNT	
        allocate(ini_F(num_points,totalChannels))					! change Nazmul CNT	
        allocate(z(num_points,totalChannels))
        allocate(dqp(num_points,totalChannels))
        allocate(dqc(num_points,totalChannels))
        allocate(dap(num_points,totalChannels))
        allocate(dac(num_points,totalChannels))
        allocate(depth(num_points))
        allocate(sk(num_points,totalChannels))
        allocate(co(num_points))
        allocate(dx(num_points,totalChannels))
        allocate(volRemain(num_points-1,totalChannels))
        allocate(froud(num_points))
        allocate(courant(num_points-1))
        allocate(oldQ(num_points, totalChannels))
        allocate(newQ(num_points, num_time, totalChannels))			! change Nazmul CNT
        allocate(added_Q(num_points, num_time, totalChannels))	 	! change Nazmul CNT
        allocate(oldArea(num_points, totalChannels))
        allocate(newArea(num_points, totalChannels))
        allocate(oldY(num_points, totalChannels))
        allocate(newY(num_points, totalChannels))
		allocate(lateralFlow(num_points, num_time, totalChannels)) 	! change Nazmul CNT
        allocate(celerity(num_points, totalChannels))
        allocate(velocity(num_points, totalChannels)) 				! change Nazmul CNT
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

        ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		! identify minimum dx (segment length) in the network
		! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		dx = 0.
        minDx = 1e10
        do j = 1,nlinks
            ncomp= frnw_g(j,1)
            do i = 1,ncomp-1
                dx(i,j) = dx_ar_g(i,j)
            end do
            minDx = min(minDx,minval(dx(1:ncomp-1,j)))
        end do

        ! TO DO:
		! * pass initial depth as initial conditions and dont arbitrarily intialize
		z=z_ar_g
        ini_y=0.05  !* [meter]
        ini_q=0.5   !*[m^3/sec]
        oldQ = -999; oldY = -999; newQ = -999; newY = -999
        dimensionless_Cr = -999; dimensionless_Fo = -999; dimensionless_Fi = -999
        dimensionless_Di = -999; dimensionless_Fc = -999; dimensionless_D = -999

		! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		! create channel geometry lookup tables
		! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		do j = 1,nlinks
			ncomp= frnw_g(j,1)
			do i=1,ncomp
			
				leftBank(i,j)= (twcc_ar_g(i,j)-tw_ar_g(i,j))/2.0
				rightBank(i,j)= (twcc_ar_g(i,j)-tw_ar_g(i,j))/2.0 + tw_ar_g(i,j)
				skLeft(i,j)= 1.0/manncc_ar_g(i,j)
				skRight(i,j)= 1.0/manncc_ar_g(i,j)
				skMain(i,j)= 1.0/mann_ar_g(i,j)
				
				call readXsection(i,(1.0/skLeft(i,j)),(1.0/skMain(i,j)),(1.0/skRight(i,j)),&
									leftBank(i,j), rightBank(i,j),timesDepth, j,&
									z_ar_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g)
									
				oldY(i,j) = ini_y(j) + z(i,j)
									
			enddo
		enddo

        ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		! Populate lateral inflow and upper boundary time arrays
		! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		! lateral inflow time
        do n=1, nts_ql_g
            tarr_ql(n)= t0_g*60.0 + dt_ql_g*real(n-1,KIND(dt_ql_g))/60.0 !* [min]
        end do
        ! upstream boundary time
        do n=1, nts_ub_g
            tarr_ub(n)= t0_g*60.0 + dt_ub_g*real(n-1,KIND(dt_ub_g))/60.0 !* [min]
        enddo

        ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		! Define initial parameters for Diffusive Wave
		! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
		min_Q = 0.028316    ! 1 cfs = 0.028316 m3/s
		
		ini_E = 1.0
		ini_F = 0.0
		
		! Define initial discharge conditons
		ini_q_repeat = iniq

		! initialization of Q, celerity, and diffusivity
		celerity = 1.0
		diffusivity = 100.0
		do j=1,nlinks
			ncomp = frnw_g(j,1)
			do i=1, ncomp
				newQ(i,1,j) = iniq(i,j)
			end do
		end do
		
		do kkk = 1,totalTimeSteps-1, repeatInterval

			! first time of kkk-th repeatInterval (minutes)
			ini_time = real(kkk-1)*dtini/60.+t0*60.
			
			! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			! FLOW CALCULATION
			! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			do j = 1, nlinks

				ncomp = frnw_g(j,1)

				do timestep=1, ntim
				
					! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					! applying the upstream boundary conditions in the matrix newQ
					! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					
					if (frnw_g(j,3) == 0) then ! reach is a headwater - use upstream boundary condition data
						varr_ub = ubcd_g(:,j)
						! interpolate upper boundary data to simulation time interval
						newQ(1,timestep,j)= intp_y(nts_ub_g, tarr_ub, varr_ub, ini_time+dtini/60.*(timestep-1)) !* tarr_ub in min.
					end if
					
					if (frnw_g(j,3).gt.0) then  ! reach boundary originates from a junction
					
						newQ(1,timestep,j) = 0.
						do k=1, frnw_g(j,3)
						
							linknb= frnw_g(j,3+k)
							nodenb= frnw_g(linknb,1)
							
							! sum of junction inflows
							newQ(1,timestep,j) = newQ(1,timestep,j) + newQ(nodenb,timestep,linknb)
		
						end do

					end if
				  
					! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					! populate lateral inflow array
					! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					
					lateralFlow(:,:,j) = 0
					do i=1,ncomp-1
						do n=1,nts_ql_g
							varr_ql(n)= qlat_g(n,i,j)
						enddo
						lateralFlow(i,timestep,j)= intp_y(nts_ql_g, tarr_ql, varr_ql, ini_time+dtini/60.*(timestep-1)) ! tarr_ql in minutes
					end do
					
				end do

				! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
				! Diffusive wave forward sweep to calculate flow at t+1
				! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					
				qp(:,:,j) = newQ(:,:,j)
				call diffusive_CNT(j,ntim,repeatInterval)

				newQ(:,:,j) = qp(:,:,j)

			end do

			! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			! Write flow results to q_ev_g array
			! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			
			do timestep = 1,repeatInterval,saveFrequency
				ts_ev = ((ini_time*60.0)/dtini) + (timestep-1)
				do j=1,nlinks
					ncomp = frnw_g(j,1)
					do i=1, ncomp
						q_ev_g(ts_ev, i, j)= newQ(i,timestep,j)
					end do
				end do
			end do
			
			! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			! DEPTH, CELERITY, and DIFFUSIVITY CALCULATION
			! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			added_Q = 0.

			do j=nlinks,1,-1

				ncomp = frnw_g(j,1)
				
				! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
				! applying the downstream boundary conditons in matrix newY
				! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
				if (frnw_g(j,2) < 0.0) then    ! use normal depth water depth calculation
				
					q_sk_multi=1.0
					slope = (z(ncomp-1,j)-z(ncomp,j))/dx(ncomp-1,j)
					if (slope .le. 0.0001) slope = 0.0001

					! calculate normal depth 
					call normal_crit_y(ncomp, j, q_sk_multi, slope, newQ(ncomp,repeatInterval,j), oldY(ncomp,j), temp,  temp, temp)
					
					
				elseif (frnw_g(j,2) .ge. 0.0) then    ! water level is calculated from the downstream river reach
					linknb=frnw_g(j,2) ! reach downstream of j
					newY(ncomp,j)= newY(1,linknb)   ! taking the WL from the d/s reach
				end if

				call mesh_diffusive_backward(j,ntim,repeatInterval)
				
			end do
			
			! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			! Write depth results to q_ev_g array
			! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			do timestep = 1,repeatInterval,saveFrequency
				ts_ev = ((ini_time*60.0)/dtini) + (timestep-1)
				do j=1,nlinks
					ncomp = frnw_g(j,1)
					do i=1, ncomp
						elv_ev_g(ts_ev, i, j)= newY(i,j) ! this will just repeat the deth value for all steps in the repeatInterval
					end do
				end do
			end do

			oldY = newY

		end do ! end kkk loop
	end subroutine diffnw
	
! ###########################################################################################
! ###########################################################################################
! ###########################################################################################

	subroutine diffusive_CNT(j,ntim,repeatInterval)

		use constants_module
		use arrays_module 
		use matrix_module 
		use var_module	  
		use arrays_section_module 
		use xsec_attribute_module 
		use subtools

		implicit none

		integer, intent(in) :: j, ntim, repeatInterval

		real,allocatable :: E_cnt(:,:), F_cnt(:,:)

		integer :: i, n, kk

		real :: hi, gi, ki, pj, qj, rj, pj_p, qj_p, rj_p, sj_p, mi, ni, qp_ghost, r_interpol_time, qp_ghost_1


	!!++++++++++++++++++++ Diffusive wave Forward sweep starts +++++++++++++++++++!!

		ncomp = frnw_g(j,1)

		allocate(E_cnt(ncomp,ntim))
		allocate(F_cnt(ncomp,ntim))

		qp_ghost_1 = qp(1,ntim,j)

		E_cnt(1:ncomp,1) = ini_E(1:ncomp,j)
		F_cnt(1:ncomp,1) = ini_F(1:ncomp,j)

		do i = 2,ncomp

			do n = 2, ntim

				hi = dx(i-1,j) / (dtini * celerity(i-1,j))
				gi = diffusivity(i-1,j) * dx(i-1,j) / ((celerity(i-1,j)**3.0)*(dtini**2.0))
				ki = (gi ** 2.0)/(hi ** 2.0)
				mi = diffusivity(i-1,j) / (celerity(i-1,j)**2.0) * dx(i-1,j) / dtini
				ni = dx(i-1,j)

				pj = - hi / 4.0 - gi / 2.0 - 2.0 * ki
				qj = 1.0 + gi + 4.0 * ki
				rj = hi / 4.0 - gi / 2.0 - 2.0 * ki

				pj_p = hi / 4.0 + gi / 2.0 - 2.0 * ki
				qj_p = 1.0 - gi + 4.0 * ki
				rj_p = -hi / 4.0 + gi / 2.0 - 2.0 * ki

				sj_p = pj_p * qp(i-1,n-1,j) + qj_p * qp(i-1,n,j) + rj_p * qp(i-1,n+1,j) &
					+ mi / 2.0 * (lateralFlow(i-1,n+1,j) - lateralFlow(i-1,n-1,j) ) &
					+ ni * lateralFlow(i-1,n,j)

				if (n .eq. ntim) then
					! applying ghost node !&
					sj_p = pj_p * qp(i-1,n-1,j) + qj_p * qp(i-1,n,j) + rj_p * qp_ghost_1 &
					+ mi / 2.0 * (lateralFlow(i-1,n,j) - lateralFlow(i-1,n-1,j) ) &
					+ ni * lateralFlow(i-1,n,j)
				end if

				E_cnt(i,n) = -1.0 * rj / (pj * E_cnt(i,n-1) + qj)
				F_cnt(i,n) = ( sj_p - pj * F_cnt(i,n-1) ) / ( pj * E_cnt(i,n-1) + qj )

			end do

			qp_ghost = qp(i-1,ntim,j)+lateralFlow(i-1,ntim,j)*dx(i-1,j)

			qp_ghost_1 = qp_ghost

			! boundary at t=ntim, calculated from the ghost node at t=ntim+1
			qp(i,ntim,j) = E_cnt(i,ntim) * qp_ghost + F_cnt(i,ntim)

			do n = ntim-1,1,-1
				qp(i,n,j) = E_cnt(i,n) * qp(i,n+1,j) + F_cnt(i,n)

				if (qp(i,n,j) .lt. min_Q) then
					added_Q(i,n,j) = min_Q - qp(i,n,j)
					qp(i,n,j) = max(qp(i,n,j),min_Q)
				end if

			end do

		end do

		! replacing with the initial value at all nodes
		qp(1:ncomp,1,j) = ini_q_repeat(1:ncomp,j)

		! taking the new initial value for the next cycle
		ini_E(1:ncomp,j) = E_cnt(1:ncomp,repeatInterval+1)
		ini_F(1:ncomp,j) = F_cnt(1:ncomp,repeatInterval+1)
		ini_q_repeat(1:ncomp,j) = qp(1:ncomp,repeatInterval+1,j)

		deallocate (E_cnt, F_cnt)

	end subroutine diffusive_CNT

! ###########################################################################################
! ###########################################################################################
! ###########################################################################################

	subroutine mesh_diffusive_backward(j,ntim,repeatInterval)

		use constants_module
		use arrays_module
		use matrix_module
		use var_module
		use arrays_section_module
		use xsec_attribute_module
		use subtools

		implicit none

		integer, intent(in) :: j,ntim,repeatInterval

		real :: a1, a2, a3, a4, b1, b2, b3, b4, dd1, dd2, dd3, dd4, h1, h2, h3, h4, xt
		real :: qy, qxy, qxxy, qxxxy, ppi, qqi, rri, ssi, sxi, mannings, Sb, width, slope

		real :: cour, cour2, q_sk_multi, sfi, r_interpol_time, r_interpo_nn, temp, dkdh
		real :: D_lim1, D_lim2, y_norm, y_crit, area_n, area_c, chnWidth, vel,y_norm1

		real :: y_norm_ds, y_crit_ds, S_ncomp, frds, area_0, width_0, hydR_0, errorY, currentQ, stg1, stg2
		integer :: tableLength, jj, iii

		real :: elevTable_1(nel),areaTable_1(nel),rediTable_1(nel),convTable_1(nel),topwTable_1(nel),currentSquaredDepth_1(nel)
		real :: dKdATable_1(nel)
		real :: pereTable_1(nel),depthYi,tempDepthi_1,tempCo_1,temp_q_sk_multi_1,tempY_1,tempArea_1,tempRadi_1,tempbo_1,ffy
		real :: ffy_1, ffy_2, ffy_3, tempCo_2, tempCo_3, tempsfi_2, tempsfi_3
		real :: ffprime,tempDepthi_1_new,tempsfi_1,toll, dkda, tempPere_1, tempdKdA_1, tempY_2, tempY_3, temp_v_1, tempDepthi_2
		real :: skkkTable_1(nel), tempsk_1
		real :: usFroud, dsFroud, dUdt, eHds, y_alt, area_alt, y_cnj, area_cnj, y_crit_test, y_norm_test
		integer :: depthCalOk(ncomp), wlCalcMethod

		integer :: i, pp

	!!++++++++++++++++++++ Diffusive wave Backward sweep starts +++++++++++++++++++!!

			D_lim1 = -10.
			D_lim2 = -15.
			wlCalcMethod = 3

			! wlCalcMethod = 1 : Mid-point bisection method
			! wlCalcMethod = 2 : simple iterative method method with contribution from dU/dX
			! wlCalcMethod = 3 : Newton Raphson method
			! wlCalcMethod = 6 : Only Normal depth

			ncomp = frnw_g(j,1)

			depthCalOk(ncomp) = 1

			do i=ncomp,1,-1

				currentQ = qp(i,repeatInterval+1,j)
				call calc_q_sk_multi(i,j,currentQ,q_sk_multi)

				! Calculating : read all attributes from tab file
				elevTable = xsec_tab(1,:,i,j)
				convTable = xsec_tab(5,:,i,j)
				areaTable = xsec_tab(2,:,i,j)
				pereTable = xsec_tab(3,:,i,j)
				topwTable = xsec_tab(6,:,i,j)
				skkkTable = xsec_tab(11,:,i,j)
				
				! interpolate the cross section attributes based on water elevation
				xt=newY(i,j)

				currentSquareDepth=(elevTable-z(i,j))**2.

				call r_interpol(currentSquareDepth,convTable,nel,(newY(i,j)-z(i,j))**2.0,co(i))
				if (co(i) .eq. -9999) then
					print*, 'At j = ',j,', i = ',i, 'interpolation of conveyence was not possible, wl', &
					newY(i,j), 'z',z(i,j),'previous wl',newY(i+1,j), 'previous z',z(i+1,j)
					stop
				end if
				co(i) =q_sk_multi * co(i)

				call r_interpol(elevTable,areaTable,nel,xt,newArea(i,j))
				if (newArea(i,j) .eq. -9999) then
					print*, 'At j = ',j,', i = ',i, 'interpolation of newArea was not possible'
					stop
				end if
				call r_interpol(elevTable,pereTable,nel,xt,pere(i,j))
				call r_interpol(elevTable,topwTable,nel,xt,bo(i,j))
				call r_interpol(elevTable,skkkTable,nel,xt,sk(i,j))

				sfi = qp(i,repeatInterval+1,j) * abs(qp(i,repeatInterval+1,j)) / ( co(i)** 2.0 )
				if (abs(qp(i,repeatInterval+1,j)) .lt. min_Q) then
					sfi = min_Q ** 2.0 * sign(1.0,qp(i,repeatInterval+1,j)) / ( co(i)** 2.0 )
					! at some head water basin, the Q boundary becomes 0.0 m3/s at some time. If Q = 0, then
					! sfi = 0. This leads to diffusivity = NaN. To avoid NaN diffusivity, those sfi is calculated using the min_Q
					! for those cases.
				end if

				chnWidth = rightBank(i,j)-leftBank(i,j)
				chnWidth = min(chnWidth,bo(i,j))

				if (depthCalOk(i) .eq. 1) then
					celerity2(i)=5.0 / 3.0 * abs(sfi) ** 0.3 * abs(qp(i,repeatInterval+1,j)) ** &
						0.4 / bo(i,j) ** 0.4 / (1./(sk(i,j)*q_sk_multi)) ** 0.6
					diffusivity2(i) = abs(qp(i,repeatInterval+1,j)) / 2.0 / bo(i,j) / abs(sfi)
					vel = qp(i,repeatInterval+1,j)/newArea(i,j)

					velocity(i,j) = vel

					! Applying the upper limit of celerity
					if (celerity2(i) .gt. 3.0*vel) celerity2(i) = vel*3.0
				else
					if (qp(i,repeatInterval+1,j) .lt. 1) then
						celerity2(i)=0.5
					else
						celerity2(i)=1.0
					end if

					diffusivity2(i)=diffusivity(i,j)
				end if

				! Calculating Depth at the upstream node
				if (i .gt. 1) then

					slope = (z(i-1,j)-z(i,j))/dx(i-1,j)
					depthYi = newY(i,j) - z(i,j)

					tempDepthi_1 = oldY(i-1,j)-z(i-1,j)

					elevTable_1 = xsec_tab(1,:,i-1,j)
					areaTable_1 = xsec_tab(2,:,i-1,j)
					pereTable_1 = xsec_tab(3,:,i-1,j)
					rediTable_1 = xsec_tab(4,:,i-1,j)
					convTable_1 = xsec_tab(5,:,i-1,j)
					topwTable_1 = xsec_tab(6,:,i-1,j)
					dKdATable_1 = xsec_tab(9,:,i-1,j)

					currentSquaredDepth_1=(elevTable_1-z(i-1,j))**2.

					toll = 1.0
					iii = 0

					dsFroud = abs(qp(i,repeatInterval+1,j))/sqrt(grav*newArea(i,j)**3.0/bo(i,j))

					! Applying mid point bisection
					if (wlCalcMethod .eq. 1) then
						tempY_1 = elevTable_1(2)
						tempY_2 = depthYi * 3. + z(i-1,j)
						tempY_3 = (tempY_1 + tempY_2) / 2.
						do while ( abs(toll) .gt. 0.001)
							iii = iii +1

							call r_interpol(currentSquaredDepth_1,convTable_1,nel,(tempY_1-z(i-1,j))**2.0,tempCo_1)
							call r_interpol(currentSquaredDepth_1,convTable_1,nel,(tempY_2-z(i-1,j))**2.0,tempCo_2)
							call r_interpol(currentSquaredDepth_1,convTable_1,nel,(tempY_3-z(i-1,j))**2.0,tempCo_3)

							call calc_q_sk_multi(i-1,j,qp(i-1,repeatInterval+1,j),temp_q_sk_multi_1)
							tempCo_1 = tempCo_1 * temp_q_sk_multi_1
							tempCo_2 = tempCo_2 * temp_q_sk_multi_1
							tempCo_3 = tempCo_3 * temp_q_sk_multi_1

							tempsfi_1 = qp(i-1,repeatInterval+1,j) * abs(qp(i-1,repeatInterval+1,j)) / ( tempCo_1** 2.0 )
							tempsfi_2 = qp(i-1,repeatInterval+1,j) * abs(qp(i-1,repeatInterval+1,j)) / ( tempCo_2** 2.0 )
							tempsfi_3 = qp(i-1,repeatInterval+1,j) * abs(qp(i-1,repeatInterval+1,j)) / ( tempCo_3** 2.0 )

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

					! Applying simple iteration method with inclusion of dU/dX
					if (wlCalcMethod .eq. 2) then
						vel = qp(i,repeatInterval+1,j)/newArea(i,j)
						toll = 1.0

						do while ( abs(toll) .gt. 0.001)
							iii = iii +1
							tempY_1 = tempDepthi_1 + z(i-1,j)

							call r_interpol(currentSquaredDepth_1,convTable_1,nel,(tempDepthi_1)**2.0,tempCo_1)

							call calc_q_sk_multi(i-1,j,qp(i-1,repeatInterval+1,j),temp_q_sk_multi_1)
							tempCo_1 = tempCo_1 * temp_q_sk_multi_1

							call r_interpol(elevTable_1,areaTable_1,nel,tempY_1,tempArea_1)

							temp_v_1 = qp(i-1,repeatInterval+1,j) / tempArea_1

							tempsfi_1 = qp(i-1,repeatInterval+1,j) * abs(qp(i-1,repeatInterval+1,j)) / ( tempCo_1** 2.0 )!

							tempDepthi_2 = depthYi - dx(i-1,j)*slope + 0.5*dx(i-1,j)*(sfi + tempsfi_1)&
								+0.5*(vel+temp_v_1)*(vel-temp_v_1)/grav

							toll = tempDepthi_2 - tempDepthi_1

							tempDepthi_1 = tempDepthi_2

						end do
						depthCalOk(i-1) = 1

					end if

					! Applying NewtonRaphson method corrected
					if (wlCalcMethod .eq. 3) then
						vel = qp(i,repeatInterval+1,j)/newArea(i,j)
						do while ( abs(toll) .gt. 0.001)
							iii = iii +1
							tempY_1 = tempDepthi_1 + z(i-1,j)

							call r_interpol(currentSquaredDepth_1,convTable_1,nel,(tempDepthi_1)**2.0,tempCo_1)
							call calc_q_sk_multi(i-1,j,qp(i-1,repeatInterval+1,j),temp_q_sk_multi_1)
							tempCo_1 = tempCo_1 * temp_q_sk_multi_1

							call r_interpol(elevTable_1,areaTable_1,nel,tempY_1,tempArea_1)
							call r_interpol(elevTable_1,pereTable_1,nel,tempY_1,tempPere_1)
							call r_interpol(elevTable_1,rediTable_1,nel,tempY_1,tempRadi_1)
							call r_interpol(elevTable_1,topwTable_1,nel,tempY_1,tempbo_1)
							call r_interpol(elevTable_1,dKdATable_1,nel,tempY_1,tempdKdA_1)!    Change 20210520

							tempsfi_1 = qp(i-1,repeatInterval+1,j) * abs(qp(i-1,repeatInterval+1,j)) / ( tempCo_1** 2.0 )!

							temp_v_1 = qp(i-1,repeatInterval+1,j)/tempArea_1

							ffy = tempDepthi_1- depthYi + dx(i-1,j)*slope - 0.5*dx(i-1,j)*(sfi + tempsfi_1)

							dkda = tempdKdA_1

							ffprime = 1 + dx(i-1,j) * tempbo_1 *  qp(i-1,repeatInterval+1,j) * &
								abs(qp(i-1,repeatInterval+1,j)) / (tempCo_1 ** 3.0) * dkda

							if (ffprime .eq. 0.) then
								! print*, 'ffprime = 0'
								! print*, j, i, qp(i-1,repeatInterval+1,j), tempCo_1
								! pause
							end if

							tempDepthi_1_new = tempDepthi_1 - ffy / ffprime

							tempDepthi_1_new = max(tempDepthi_1_new,0.005)

							toll = abs(tempDepthi_1_new - tempDepthi_1)

							if(iii .gt. 30)then
								print*, 'Warning: Depth iteration reached maximum trial at j=', j, 'i=', i-1 , 'and',i, &
								'depths are',tempDepthi_1, tempDepthi_1_new, 'slope=', slope, 'tempbo_1=',tempbo_1, 'dx=',&
								dx(i-1,j), 'Q=', qp(i-1,repeatInterval+1,j)
								print*, 'depth at d/s',depthYi
								depthCalOk(i-1) = 0
								EXIT
							endif
							tempDepthi_1 = tempDepthi_1_new
							depthCalOk(i-1) = 1
						end do

					end if      ! calculation of Wl based on d/s wl ends


					usFroud = abs(qp(i-1,repeatInterval+1,j))/sqrt(grav*tempArea_1**3.0/tempbo_1)

					newY(i-1,j) = tempDepthi_1 + z(i-1,j)

					! calculating the normal depth at i=i-1 using the slope between i and i-1
					! Calculating the normal depth and critical depth at the river reach upstream as an output
					! very small slope is creating an issue. Applying a min limit of bed slope.
					! for Normal depth calculation, applying the absolute value of the Q
					call normal_crit_y(i-1, j, q_sk_multi, max(slope,0.0001), max(abs(qp(i-1,repeatInterval+1,j)),min_Q),&
						y_norm, y_crit, area_norm, area_crit)

					! to avoid sudden jump in water level calculation
					! if the WL jumps abruptly, the tempsfi_1 would be very small compared to sfi.
					! we are replacing the water level using normal water level
					if (sfi / tempsfi_1 .gt. 1e5)  newY(i-1,j) = y_norm

					! limiting very high value of WL
					if (newY(i-1,j) .gt. 10.0**5.) newY(i-1,j) = 10.0**5.

					if (newY(i-1,j)-z(i-1,j) .le. 0.) then
						! print*, 'depth is negative at j= ', j,'i=',i-1,'newY=',(newY(jj,j),jj=1,ncomp)
						! print*, 'dimensionless_D',(dimensionless_D(jj,j),jj=1,ncomp)
						! print*, 'newQ',(qp(jj,repeatInterval+1,j),jj=1,ncomp)
						! print*, 'Bed',(z(jj,j),jj=1,ncomp)
						! print*, 'dx',(dx(jj,j),jj=1,ncomp-1)
						! pause 777
					end if

				end if      ! end of if (i .gt. 1) || end of WL calculation at j reach

				! Book-keeping: Counting the number as for how many time steps the routing method is unchanged
				routingNotChanged(i-1,j) = routingNotChanged(i-1,j) + 1

			end do

			celerity(1:ncomp,j) =  sum(celerity2(1:ncomp)) / ncomp  ! change 20210524

			do i = 1, ncomp
				celerity(i,j) = max(celerity(i,j),0.5)   !!! Applying Celerity lower limit
			end do

			diffusivity(1:ncomp,j)=sum(diffusivity2(1:ncomp)) / ncomp

			do i = 1, ncomp
				if (diffusivity(i,j) .gt. maxDiffuLm) diffusivity(i,j) = maxDiffuLm !!! Applying diffusivity upper limit
				if (diffusivity(i,j) .lt. minDiffuLm) diffusivity(i,j) = minDiffuLm !!! Applying diffusivity lower limit
			end do

	!!++++++++++++++++++++ Diffusive wave Backward sweep ends +++++++++++++++++++!!

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
	
endmodule diffusive





























