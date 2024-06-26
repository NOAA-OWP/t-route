!-------------------------------------------------------------------------------
module diffusive

  IMPLICIT NONE
  
!-----------------------------------------------------------------------------
! Description:
!   Numerically solve diffusive wave PDEs using Crank-Nicolson and Hermite
!   Interpolation. 
!
! Current Code Owner: NOAA-OWP, Inland Hydraulics Team
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
  
  ! Module constants
  double precision, parameter :: grav = 9.81
  double precision, parameter :: TOLERANCE = 1e-8
  
  ! Module variables
  integer :: nlinks
  integer :: mxncomp
  integer :: maxTableLength
  integer :: nel
  integer :: newtonRaphson
  integer :: nmstem_rch
  integer :: nts_da
  integer, dimension(:),   allocatable :: mstem_frj    
  integer, dimension(:),   allocatable :: usgs_da_reach  
  integer, dimension(:,:), allocatable :: frnw_g
  double precision :: dtini, dxini, cfl, minDx, maxCelerity,  theta
  double precision :: C_llm, D_llm, D_ulm, DD_ulm, DD_llm, q_llm, so_llm
  double precision :: dt_qtrib
  double precision :: dmyi, dmyj
  double precision :: dtini_min
  double precision, dimension(:),       allocatable :: eei, ffi, exi, fxi, qcx, diffusivity2, celerity2
  double precision, dimension(:),       allocatable :: ini_y, ini_q
  double precision, dimension(:),       allocatable :: elevTable, areaTable, skkkTable
  double precision, dimension(:),       allocatable :: pereTable, rediTable
  double precision, dimension(:),       allocatable :: convTable, topwTable
  double precision, dimension(:),       allocatable :: currentSquareDepth
  double precision, dimension(:),       allocatable :: tarr_qtrib, varr_qtrib
  double precision, dimension(:),       allocatable :: tarr_da, varr_da
  double precision, dimension(:),       allocatable :: area, depth, co, froud, courant
  double precision, dimension(:,:),     allocatable :: bo, dx
  double precision, dimension(:,:),     allocatable :: areap, qp, z, sk
  double precision, dimension(:,:),     allocatable :: celerity, diffusivity, qpx  
  double precision, dimension(:,:),     allocatable :: pere, oldQ, newQ, oldArea, newArea, oldY, newY  
  double precision, dimension(:,:),     allocatable :: volRemain  
  double precision, dimension(:,:),     allocatable :: lateralFlow
  double precision, dimension(:,:),     allocatable :: dimensionless_Cr, dimensionless_Fo, dimensionless_Fi
  double precision, dimension(:,:),     allocatable :: dimensionless_Di, dimensionless_Fc, dimensionless_D  
  double precision, dimension(:,:),     allocatable :: qtrib 
  double precision, dimension(:,:),     allocatable :: usgs_da

  double precision, dimension(:,:,:,:), allocatable :: xsec_tab
  
contains

  subroutine diffnw(timestep_ar_g, nts_ql_g, nts_ub_g, nts_db_g, ntss_ev_g, nts_qtrib_g, nts_da_g,   &
                    mxncomp_g, nrch_g, dx_ar_g, iniq, frnw_col, frnw_ar_g, qlat_g, dbcd_g, qtrib_g,  &                       &
                    paradim, para_ar_g, usgs_da_g, usgs_da_reach_g, q_ev_g, elv_ev_g, depth_ev_g)                                     
                    

    IMPLICIT NONE
          
  !-----------------------------------------------------------------------------
  ! Description:
  !   Compute diffusive routing on National Water Model channel network domain.
  !
  ! Method:
  !   A Crank-Nicolson solution of the diffusive wave equations is solved with
  !   adaptive timestepping to maintain numerical stability. Major operations are
  !   as follows:
  !
  !     1. Initialize domain flow and depth. Depth is computed from initial flow
  !     2. For each timestep....
  !       2.1. Compute domain flow, upstream-to-downstream
  !       2.2. Compute domain depth, downstream-to-upstream
  !       2.3. Determine time step duration needed for stability
  !       2.4. repeat until until final time is reached
  !     3. Record results in output arrays at user-specified time intervals
  !
  ! Current Code Owner: NOAA-OWP, Inland Hydraulics Team
  !
  ! Development Team:
  !  - Dong Ha Kim
  !  - Nazmul Azim Beg
  !  - Ehab Meselhe
  !  - Adam N. Wlostowski
  !  - James Halgren
  !  - Jacob Hreha
  !
  ! Code Description:
  !   Language: Fortran 90.
  !   This code is written to JULES coding standards v1.
  !-----------------------------------------------------------------------------

  ! Subroutine arguments
    integer, intent(in) :: mxncomp_g 
    integer, intent(in) :: nrch_g
    integer, intent(in) :: nts_ql_g
    integer, intent(in) :: nts_ub_g
    integer, intent(in) :: nts_db_g
    integer, intent(in) :: ntss_ev_g
    integer, intent(in) :: nts_qtrib_g
    integer, intent(in) :: nts_da_g
    integer, intent(in) :: frnw_col
    integer, intent(in) :: paradim
    integer, dimension(nrch_g), intent(in) :: usgs_da_reach_g
    integer, dimension(nrch_g, frnw_col),  intent(in) :: frnw_ar_g
    double precision, dimension(paradim ), intent(in) :: para_ar_g
    double precision, dimension(:)       , intent(in) :: timestep_ar_g(10)
    double precision, dimension(nts_db_g), intent(in) :: dbcd_g
    double precision, dimension(mxncomp_g,   nrch_g),            intent(in) :: dx_ar_g  
    double precision, dimension(mxncomp_g,   nrch_g),            intent(in) :: iniq
    double precision, dimension(nts_qtrib_g, nrch_g),            intent(in) :: qtrib_g 
    double precision, dimension(nts_da_g,    nrch_g),            intent(in) :: usgs_da_g
    double precision, dimension(nts_ql_g,  mxncomp_g, nrch_g),   intent(in ) :: qlat_g
    double precision, dimension(ntss_ev_g, mxncomp_g, nrch_g),   intent(out) :: q_ev_g
    double precision, dimension(ntss_ev_g, mxncomp_g, nrch_g),   intent(out) :: elv_ev_g
    double precision, dimension(ntss_ev_g, mxncomp_g, nrch_g),   intent(out) :: depth_ev_g

  ! Local variables    
    integer :: ncomp
    integer :: i
    integer :: j
    integer :: k
    integer :: n
    integer :: timestep
    integer :: kkk
    integer :: ts_ev
    integer :: jm
    integer :: usrchj
    integer :: linknb
    integer :: xcolID
    integer :: ycolID
    integer :: iel
    integer :: dsbc_option
    integer :: ts, cwrow 
    integer :: ri, rj, oi, oj
    integer :: nlnk, lnk
    integer :: nusrch
    integer, dimension(:), allocatable   :: dmy_frj
    double precision :: x
    double precision :: saveInterval, width
    double precision :: maxCourant
    double precision :: dtini_given, dtini_divisor
    double precision :: timesDepth
    double precision :: t
    double precision :: tfin
    double precision :: t0
    double precision :: q_sk_multi
    double precision :: maxCelDx
    double precision :: slope
    double precision :: y_norm
    double precision :: temp
    double precision :: dt_ql
    double precision :: dt_ub
    double precision :: dt_db
    double precision :: dt_da
    double precision :: wdepth
    double precision :: q_usrch
    double precision :: tf0
    double precision :: convey
    double precision :: equiv_one

    double precision :: mindepth_nstab
    double precision, dimension(:), allocatable :: tarr_ql
    double precision, dimension(:), allocatable :: varr_ql
    double precision, dimension(:), allocatable :: tarr_ub
    double precision, dimension(:), allocatable :: varr_ub
    double precision, dimension(:), allocatable :: tarr_db
    double precision, dimension(:), allocatable :: varr_db
    double precision, dimension(:), allocatable :: dbcd 

    double precision, dimension(:,:,:), allocatable :: temp_q_ev_g
    double precision, dimension(:,:,:), allocatable :: temp_elv_ev_g    

  !-----------------------------------------------------------------------------
  ! Time domain parameters
    dtini         = timestep_ar_g(1)    ! initial timestep duration [sec]
    t0            = timestep_ar_g(2)    ! simulation start time [hr]
    tfin          = timestep_ar_g(3)    ! simulation end time [hr]
    saveInterval  = timestep_ar_g(4)    ! output recording interval [sec]
    dt_ql         = timestep_ar_g(5)    ! lateral inflow data time step [sec]
    dt_ub         = timestep_ar_g(6)    ! upstream boundary time step [sec]
    dt_db         = timestep_ar_g(7)    ! downstream boundary time step [sec]
    dt_qtrib      = timestep_ar_g(8)    ! tributary data time step [sec]
    dt_da         = timestep_ar_g(9)    ! data assimilation data time step [sec]
    dtini_given   = dtini               ! preserve user-input timestep duration
    dtini_divisor = timestep_ar_g(10)   ! numeric value that divide dtini in the next line
    dtini_min     = dtini/dtini_divisor ! minimum increment in internal time step 

  !-----------------------------------------------------------------------------
  ! miscellaneous parameters
    timesDepth = 4.0 ! water depth multiplier used in readXsection
    nel        = 501 ! number of sub intervals in look-up tables
    nts_da     = nts_da_g ! make DA time steps global

  !-----------------------------------------------------------------------------
  ! Network mapping (frnw_g) array size parameters
    mxncomp = mxncomp_g ! maximum number of nodes in a single reach
    nlinks  = nrch_g    ! number of reaches in the network

  !-----------------------------------------------------------------------------
  ! Some essential parameters for Diffusive Wave
    cfl         = para_ar_g(1)  ! maximum Courant number (default: 0.95)
    C_llm       = para_ar_g(2)  ! lower limit of celerity (default: 0.5)
    D_llm       = para_ar_g(3)  ! lower limit of diffusivity (default: 50)
    D_ulm       = para_ar_g(4)  ! upper limit of diffusivity (default: 1000)
    DD_llm      = para_ar_g(5)  ! lower limit of dimensionless diffusivity (default -15)
    DD_ulm      = para_ar_g(6)  ! upper limit of dimensionless diffusivity (default: -10.0)
    q_llm       = para_ar_g(8)  ! lower limit of discharge (default: 0.02831 cms)
    so_llm      = para_ar_g(9)  ! lower limit of channel bed slope (default: 0.0001)
    theta       = para_ar_g(10) ! weight for computing 2nd derivative: 
                                ! 0: explicit, 1: implicit (default: 1.0)
    dsbc_option = para_ar_g(11) ! downstream water depth boundary condition option 1:given water depth data, 2:normal depth   
    mindepth_nstab = 0.1        ! minimum water depth for numerical stability when estimated/observed downstream boundary depth data is used.


  !-----------------------------------------------------------------------------
  ! consider moving variable allocation to a separate module
    allocate(area(mxncomp))
    allocate(bo(mxncomp,nlinks))
    allocate(pere(mxncomp,nlinks))
    allocate(areap(mxncomp,nlinks))
    allocate(qp(mxncomp,nlinks))
    allocate(z(mxncomp,nlinks))
    allocate(depth(mxncomp))
    allocate(sk(mxncomp,nlinks))
    allocate(co(mxncomp))
    allocate(dx(mxncomp,nlinks))
    allocate(volRemain(mxncomp-1,nlinks))
    allocate(froud(mxncomp))
    allocate(courant(mxncomp-1))
    allocate(oldQ(mxncomp, nlinks))
    allocate(newQ(mxncomp, nlinks))
    allocate(oldArea(mxncomp, nlinks))
    allocate(newArea(mxncomp, nlinks))
    allocate(oldY(mxncomp, nlinks))
    allocate(newY(mxncomp, nlinks))
    allocate(lateralFlow(mxncomp,nlinks))
    allocate(celerity(mxncomp, nlinks))
    allocate(diffusivity(mxncomp, nlinks))
    allocate(celerity2(mxncomp))
    allocate(diffusivity2(mxncomp))
    allocate(eei(mxncomp))
    allocate(ffi(mxncomp))
    allocate(exi(mxncomp))
    allocate(fxi(mxncomp))
    allocate(qpx(mxncomp, nlinks))
    allocate(qcx(mxncomp))
    allocate(dimensionless_Cr(mxncomp-1,nlinks))
    allocate(dimensionless_Fo(mxncomp-1,nlinks))
    allocate(dimensionless_Fi(mxncomp-1,nlinks))
    allocate(dimensionless_Di(mxncomp-1,nlinks))
    allocate(dimensionless_Fc(mxncomp-1,nlinks))
    allocate(dimensionless_D(mxncomp-1,nlinks))
    
    allocate(elevTable(nel))
    allocate(areaTable(nel))
    allocate(pereTable(nel))
    allocate(rediTable(nel))
    allocate(convTable(nel))
    allocate(topwTable(nel))
    allocate(skkkTable(nel))
    allocate(xsec_tab(11, nel, mxncomp, nlinks))

     allocate(currentSquareDepth(nel))
    allocate(ini_y(nlinks))
    allocate(ini_q(nlinks))

    allocate(tarr_ql(nts_ql_g+1), varr_ql(nts_ql_g+1))
    allocate(tarr_ub(nts_ub_g), varr_ub(nts_ub_g))
    allocate(tarr_qtrib(nts_qtrib_g), varr_qtrib(nts_qtrib_g))
    allocate(tarr_db(nts_db_g), varr_db(nts_db_g))
    allocate(tarr_da(nts_da))
    allocate(dmy_frj(nlinks))
    allocate(frnw_g(nlinks,frnw_col))

    allocate(usgs_da_reach(nlinks))
    allocate(usgs_da(nts_da, nlinks))
    allocate(dbcd(nts_db_g))
    allocate(temp_q_ev_g(ntss_ev_g, mxncomp_g, nrch_g), temp_elv_ev_g(ntss_ev_g, mxncomp_g, nrch_g)) 

    
  !--------------------------------------------------------------------------------------------
    frnw_g        = frnw_ar_g ! network mapping matrix

    usgs_da_reach = usgs_da_reach_g ! contains indices of reaches where usgs data are available
    usgs_da       = usgs_da_g       ! contains usgs data at a related reach 
        
  !-----------------------------------------------------------------------------
  ! variable initializations
    x                   = 0.0
    newQ                = -999
    newY                = -999
    t                   = t0*60.0     ! [min]
    q_sk_multi          = 1.0
    oldQ                = iniq
    newQ                = oldQ
    qp                  = oldQ
    dimensionless_Cr    = -999
    dimensionless_Fo    = -999
    dimensionless_Fi    = -999
    dimensionless_Di    = -999
    dimensionless_Fc    = -999
    dimensionless_D     = -999
    dbcd                = dbcd_g   
    q_ev_g              = 0.0
    elv_ev_g            = 0.0
    depth_ev_g          = 0.0
    temp_q_ev_g         = 0.0
    temp_elv_ev_g       = 0.0    
  !-----------------------------------------------------------------------------
  ! Identify mainstem reaches and list their ids in an array

    ! Create a dummy array containing mainstem reaches where diffusive wave is applied.
    nmstem_rch = 0
    do j = 1, nlinks
      !if (frnw_g(j,3) >= 1) then  ! mainstem reach identification
      nusrch = frnw_g(j,3)         ! the number of upstream reaches
      if (frnw_g(j, 3+nusrch+1)==555) then ! 555 indicate diffusive reach while -555 already routed reach
        nmstem_rch          = nmstem_rch + 1
        dmy_frj(nmstem_rch) = j
      end if
    end do

    ! allocate and populate array for upstream reach ids
    allocate (mstem_frj(nmstem_rch))
    do jm = 1, nmstem_rch
      mstem_frj(jm) = dmy_frj(jm)
    end do
    deallocate(dmy_frj)

  !-----------------------------------------------------------------------------
  ! create dx array from dx_ar_g and determine minimum dx.

    dx    = 0.
    minDx = 1e10

    do jm = 1, nmstem_rch !* mainstem reach only
      j = mstem_frj(jm)
      ncomp = frnw_g(j, 1)
      do i = 1, ncomp-1
        dx(i, j) = dx_ar_g(i, j)
      end do
      minDx = min(minDx, minval(dx(1:ncomp-1, j)))
    end do
    
  


  !-----------------------------------------------------------------------------
  ! Build time arrays for lateral flow, upstream boundary, downstream boundary,
  ! and tributary flow

    ! time step series for lateral flow. 
    ! ** ql from t-route starts from the first time step. For example, when t0 and tfin = 0 and 2 hrs with dt=300 sec,
    ! ** ql values from t-route starts at 5, 10, 15, ..., 120 min.
    do n = 1, nts_ql_g
      tarr_ql(n+1) =    t0 * 60.0 + dt_ql * &
                      real(n, KIND(dt_ql)) / 60.0 ! [min]
    end do
    tarr_ql(1) = t0 * 60   ! initial time step [min]

    ! time step series for upstream boundary data
    do n = 1, nts_ub_g
      tarr_ub(n) =    t0 * 60.0 + dt_ub * &
                      real(n-1,KIND(dt_ub)) / 60.0 ! [min]
    end do

    ! time step series for tributary flow data
    ! ** qtrib from t-route starts from the initial time step . For example, when t0 and tfin = 0 and 2 hrs with dt=300 sec,
    ! ** qtrib values from t-route starts at 0, 5, 10, 15, ..., 120 min.
    do n=1, nts_qtrib_g
      tarr_qtrib(n) = t0 * 60.0 + dt_qtrib * &
                      real(n-1,KIND(dt_qtrib)) / 60.0 ! [min]
    end do
    
    ! use tailwater downstream boundary observations
    ! needed for coastal coupling
    ! **** COMING SOON ****
    ! time step series for downstream boundary data
    do n=1, nts_db_g
        tarr_db(n) = t0 * 60.0 + dt_db * &
                      real(n-1,KIND(dt_db)) / 60.0 ! [min]                 
    end do
    
    ! time step series for data assimilation for discharge at bottom node of a related reach
    do n=1, nts_da
        tarr_da(n) = t0 * 60.0 + dt_da * &
                      real(n-1,KIND(dt_da)) / 60.0 ! [min]
    end do

  !-----------------------------------------------------------------------------
  ! Initialize water surface elevation, channel area, and volume
    do jm = nmstem_rch, 1, -1
      j     = mstem_frj(jm)  ! reach index
      ncomp = frnw_g(j, 1)   ! number of nodes in reach j    
      if (frnw_g(j, 2) < 0) then 
      
        ! Initial depth at bottom node of tail water reach        
        if (dsbc_option == 1) then
        ! use tailwater downstream boundary observations
        ! needed for coastal coupling
        ! **** COMING SOON **** 
          do n = 1, nts_db_g
            varr_db(n) = dbcd(n) + z(ncomp, j) !* when dbcd is water depth [m], channel bottom elev is added.
          end do
          t              = t0 * 60.0
          oldY(ncomp, j) = intp_y(nts_db_g, tarr_db, varr_db, t)
          newY(ncomp, j) = oldY(ncomp, j)  
          if ((newY(ncomp, j) - z(ncomp, j)).lt.mindepth_nstab) then
            newY(ncomp, j) = mindepth_nstab + z(ncomp, j)
          end if          
        else if (dsbc_option == 2) then
        ! normal depth as TW boundary condition
          xcolID         = 10
          ycolID         = 1
          oldY(ncomp, j) = intp_xsec_tab(ncomp, j, nel, xcolID, ycolID, oldQ(ncomp,j)) ! normal elevation not depth
          newY(ncomp, j) = oldY(ncomp, j)  
        endif
      else
        ! Initial depth at bottom node of interior reach is equal to the depth at top node of the downstream reach
        linknb         = frnw_g(j, 2)
        newY(ncomp, j) = newY(1, linknb)        
      end if              
     
      ! compute newY(i, j) for i=1, ncomp-1 with the given newY(ncomp, j)
      ! ** At initial time, oldY(i,j) values at i < ncomp, used in subroutine rtsafe, are not defined.
      ! ** So, let's assume the depth values are all nodes equal to depth at the bottom node.
      wdepth = newY(ncomp, j) - z(ncomp, j)
      do i = 1, ncomp -1
        oldY(i,j) = wdepth + z(i, j)      
      end do
      
      call mesh_diffusive_backward(dtini_given, t0, t, tfin, saveInterval, j)

      do i = 1,ncomp      
        ! copy computed initial depths to initial depth array for first timestep
        oldY(i, j) = newY(i, j)
        
        ! Check that node elevation is not lower than bottom node.
        ! If node elevation is lower than bottom node elevation, correct.
        if (oldY(i, j) .lt. oldY(ncomp, nlinks)) oldY(i, j) = oldY(ncomp, nlinks)
      end do
      
    end do

  !-----------------------------------------------------------------------------
  ! Write tributary results to output arrays
  ! TODO: consider if this is necessary - output arrays are immediately trimmed
  ! to exclude tributary results (from MC) and pass-out only diffusive-calculated
  ! flow and depth on mainstem segments.
          
    ts_ev = 1
    do while (t <= tfin * 60.0)
      if ( (mod( (t - t0 * 60.) * 60., saveInterval) <= TOLERANCE) &
            .or. (t == tfin * 60.) ) then
        do j = 1, nlinks
          if (all(mstem_frj /= j)) then ! NOT a mainstem reach
            do n = 1, nts_qtrib_g
              varr_qtrib(n) = qtrib_g(n, j)
            end do
            q_ev_g(ts_ev, frnw_g(j, 1), j) = intp_y(nts_qtrib_g, tarr_qtrib, &
                                                      varr_qtrib, t)
            q_ev_g(ts_ev,            1, j) = q_ev_g(ts_ev, frnw_g(j, 1), j)            
          end if
        end do
        ts_ev = ts_ev + 1
      end if
      t = t + dtini / 60. !* [min]
    end do
 
  !-----------------------------------------------------------------------------
  ! Initializations and re-initializations
    qpx                     = 0.
    width                   = 100.
    maxCelerity             = 1.0
    maxCelDx                = maxCelerity / minDx
    dimensionless_Fi        = 10.1
    dimensionless_Fc        = 10.1
    dimensionless_D         = 0.1
    timestep                = 0
    ts_ev                   = 1
    t                       = t0 * 60.0

  !-----------------------------------------------------------------------------
  ! Ordered network routing computations

    do while ( t < tfin * 60.)
      timestep = timestep + 1 ! advance timestep
      !+-------------------------------------------------------------------------
      !+                             PREDICTOR
      !+
      !+    Applies only to mainstem
      !+-------------------------------------------------------------------------
      do jm = 1, nmstem_rch   ! loop over mainstem reaches [upstream-to-downstream]
        j     = mstem_frj(jm) ! reach index
        ncomp = frnw_g(j,1)   ! number of nodes in reach j

        ! Calculate the duration of this timestep (dtini)
        ! Timestep duration is selected to maintain numerical stability
        if (j == mstem_frj(1)) then 
          call calculateDT(t0, t, saveInterval, cfl, tfin, maxCelDx, dtini_given)
        end if                                        

        ! estimate lateral flow at current time t
        do i = 1, ncomp - 1
          do n = 1, nts_ql_g
            varr_ql(n+1) = qlat_g(n, i, j)
          end do
          varr_ql(1)        = qlat_g(1, i, j)
          lateralFlow(i, j) = intp_y(nts_ql_g+1, tarr_ql, varr_ql, t)
        end do

        !+++----------------------------------------------------------------
        !+ Hand over water from upstream to downstream properly according
        !+ to network connections, i.e., serial or branching.
        !+ Refer to p.52, RM1_MESH
        !+++----------------------------------------------------------------
        if (frnw_g(j, 3) > 0) then        ! number of upstream reaches, so reach j is not a headwater reach's index
          newQ(1, j) = 0.0
          do k = 1, frnw_g(j, 3)          ! loop over upstream connected reaches
            usrchj = frnw_g(j, 3 + k) 
            if (any(mstem_frj == usrchj)) then

              ! inflow from upstream mainstem reach's bottom node
              q_usrch = newQ(frnw_g(usrchj, 1), usrchj)
            else

              ! inflow from upstream tributary reach
              do n = 1, nts_qtrib_g
                varr_qtrib(n) = qtrib_g(n, usrchj)
              end do
              tf0 = t +  dtini / 60.
              q_usrch = intp_y(nts_qtrib_g, tarr_qtrib, varr_qtrib, tf0)
            end if
            ! add upstream flows to reach head
            newQ(1,j)= newQ(1,j) + q_usrch           
          end do        
        else

          ! no stream reaches at the upstream of the reach (frnw_g(j,3)==0) such as a headwater reach among tributaries
          newQ(1,j)= 0.0          
        end if

        ! Add lateral inflows to the reach head
        newQ(1, j) = newQ(1, j) + lateralFlow(1, j) * dx(1, j)

        call mesh_diffusive_forward(dtini_given, t0, t, tfin, saveInterval, j)

      end do  ! end of j loop for predictor
      
      !+-------------------------------------------------------------------------
      !+                             CORRECTOR
      !+
      !+    Applies only to mainstem
      !+-------------------------------------------------------------------------
      do jm = nmstem_rch, 1, -1 ! loop over mainstem reaches [downstream-to-upstream]
        j     = mstem_frj(jm)
        ncomp = frnw_g(j,1)
          
        !+++------------------------------------------------------------+
        !+ Downstream boundary condition for water elevation either at
        !+ a junction or TW.
        !+ Refer to p.53-1,RM1_MESH
        !+++------------------------------------------------------------+
         if (frnw_g(j, 2) >= 0.0) then          

            ! reach of index j has a reach in the downstream 
            ! set bottom node WSEL equal WSEL in top node of downstream reach
            linknb         = frnw_g(j,2)
            newY(ncomp, j) = newY(1, linknb) 
        else
        
          ! reach of index j has NO downstream connection, so it is a tailwater reach
          if (dsbc_option == 1) then  
          ! use downstream boundary data source to set bottom node WSEL value       
            newY(ncomp, j) = intp_y(nts_db_g, tarr_db, varr_db, t+dtini/60.)
            ! to prevent too small water depth at boundary
            if ((newY(ncomp, j) - z(ncomp, j)).lt.mindepth_nstab) then
              newY(ncomp, j) = mindepth_nstab + z(ncomp, j)
            end if 
            xcolID  = 1
            ycolID  = 2
            newArea(ncomp, j) = intp_xsec_tab(ncomp, j, nel, xcolID, ycolID, newY(ncomp,j)) ! area of normal elevation
          else if (dsbc_option == 2) then  
          ! Assume normal depth at tailwater downstream boundary
            xcolID  = 10
            ycolID  = 1
            newY(ncomp,j) = intp_xsec_tab(ncomp, j, nel, xcolID, ycolID, abs(newQ(ncomp,j))) ! normal elevation not depth
            xcolID  = 1
            ycolID  = 2
            newArea(ncomp,j) = intp_xsec_tab(ncomp, j, nel, xcolID, ycolID, newY(ncomp,j)) ! area of normal elevation
          endif
        end if

        ! Calculate WSEL at interior reach nodes
        call mesh_diffusive_backward(dtini_given, t0, t, tfin, saveInterval, j)
        
        ! Identify the maximum calculated celerity/dx ratio at this timestep
        ! maxCelDx is used to determine the duration of the next timestep
        if (jm == 1) then
          maxCelDx    = 0.
          do i = 1, nmstem_rch
            do kkk = 1, frnw_g(mstem_frj(i), 1)-1
              maxCelDx    = max(maxCelDx, celerity(kkk, mstem_frj(i)) / dx(kkk, mstem_frj(i)))
            end do
          end do
        endif
      end do  ! end of corrector j loop
    
      !+-------------------------------------------------------------------------
      !+                             BOOK KEEPING
      !+
      !+-------------------------------------------------------------------------
      
      ! Calculate Froude and maximum Courant number
      do jm = 1, nmstem_rch
        j     = mstem_frj(jm)
        ncomp = frnw_g(j,1)
        do i = 1, ncomp
!          froud(i) = abs(newQ(i, j)) / sqrt(grav * newArea(i, j) ** 3.0 / bo(i, j))
          if (i < ncomp) then
            courant(i) = (newQ(i, j) + newQ(i+1, j)) / (newArea(i, j) + newArea(i+1, j)) * dtini / dx(i, j)
          endif
        end do
        !if (maxCourant < maxval(courant(1:ncomp - 1))) then
        !  maxCourant = maxval(courant(1:ncomp-1))
        !end if
      end do

      ! Advance model time
      t = t + dtini/60.
      
      ! Calculate dimensionless numbers for each reach
      !do jm = 1, nmstem_rch  !* mainstem reach only
      !  j = mstem_frj(jm)
      !  call calc_dimensionless_numbers(j)
      !end do

      ! diffusive wave simulation time print
      if (mod(t,30.)==0.) then
        print*, "diffusive simulation time in minute=", t
      endif

      ! write results to output arrays
      if ( (mod((t - t0 * 60.) * 60., saveInterval) <= TOLERANCE) .or. (t == tfin * 60.)) then
        do jm = 1, nmstem_rch
          j     = mstem_frj(jm)
          ncomp = frnw_g(j, 1)
          do i = 1, ncomp
            q_ev_g  (ts_ev + 1, i, j)   = newQ(i, j)
            elv_ev_g(ts_ev + 1, i, j)   = newY(i, j)
            depth_ev_g(ts_ev + 1, i, j) = elv_ev_g(ts_ev + 1, i, j) - z(i, j)
          end do
              
          !* water elevation for tributaries flowing into the mainstem in the middle or at the upper end
          do k = 1, frnw_g(j, 3)
            usrchj = frnw_g(j, 3 + k)
            if (all(mstem_frj /= usrchj)) then
              wdepth                                         = newY(1, j) - z(1, j)
              elv_ev_g(ts_ev+1, frnw_g(usrchj, 1), usrchj)   = newY(1, j)
              depth_ev_g(ts_ev+1, frnw_g(usrchj, 1), usrchj) = wdepth
            endif
          end do
        end do
        
        ! Advance recording timestep
        ts_ev = ts_ev+1
      end if

      ! write initial conditions to output arrays
      if ( ( t == t0 + dtini / 60. ) ) then
        do jm = 1, nmstem_rch  !* mainstem reach only
          j     = mstem_frj(jm)
          ncomp = frnw_g(j, 1)
          do i = 1, ncomp
            q_ev_g  (1, i, j)   = oldQ(i, j)
            elv_ev_g(1, i, j)   = oldY(i, j)
            depth_ev_g(1, i, j) = elv_ev_g(1, i, j) - z(i, j)
          end do
              
          !* water elevation for tributaries flowing into the mainstem in the middle or at the upper end
          do k = 1, frnw_g(j, 3) !* then number of upstream reaches
            usrchj = frnw_g(j, 3 + k) !* js corresponding to upstream reaches
            if (all(mstem_frj /= usrchj)) then                  
              wdepth                                   = oldY(1, j) - z(1, j)
              elv_ev_g(1, frnw_g(usrchj,1), usrchj)    = oldY(1,j)
              depth_ev_g(1, frnw_g(usrchj, 1), usrchj) = wdepth
            end if
          end do
        end do
      end if

      ! update of Y, Q and Area vectors
      oldY    = newY
      newY    = -999
      oldQ    = newQ
      newQ    = -999
      oldArea = newArea
      newArea = -999
      pere    = -999
      
    end do  ! end of time loop
    
    deallocate(frnw_g)
    deallocate(area, bo, pere, areap, qp, z,  depth, sk, co, dx) 
    deallocate(volRemain, froud, courant, oldQ, newQ, oldArea, newArea, oldY, newY)
    deallocate(lateralFlow, celerity, diffusivity, celerity2, diffusivity2)
    deallocate(eei, ffi, exi, fxi, qpx, qcx)
    deallocate(dimensionless_Cr, dimensionless_Fo, dimensionless_Fi)
    deallocate(dimensionless_Di, dimensionless_Fc, dimensionless_D)
    deallocate(elevTable, areaTable, pereTable, rediTable, convTable, topwTable)
    deallocate(skkkTable)
    deallocate(xsec_tab)
    deallocate(currentSquareDepth, ini_y, ini_q)
    deallocate(tarr_ql, varr_ql, tarr_ub, varr_ub, tarr_qtrib, varr_qtrib, tarr_da)
    deallocate(mstem_frj)
    deallocate(usgs_da_reach, usgs_da)

  end subroutine diffnw
  
  subroutine calculateDT(initialTime, time, saveInterval, &
                         maxAllowCourantNo, tfin, max_C_dx, given_dt)
                         
    IMPLICIT NONE
  
  !-----------------------------------------------------------------------------
  ! Description:
  !   Calculate time step duration. In order to maintain numerical stability,
  !   the model timestep duration must be short enough to keep Courant numbers
  !   below a threshold level. Else, numerical instabilities will occur. 
  !
  ! Method:
  !   
  !
  ! Current Code Owner: NOAA-OWP, Inland Hydraulics Team
  !
  ! Code Description:
  !   Language: Fortran 90.
  !   This code is written to JULES coding standards v1.
  !-----------------------------------------------------------------------------

  ! Subroutine arguments
    double precision, intent(in) :: initialTime       ! [sec]
    double precision, intent(in) :: time              ! [hrs]
    double precision, intent(in) :: saveInterval      ! [sec]
    double precision, intent(in) :: tfin              ! [hrs]
    double precision, intent(in) :: given_dt          ! [sec]
    double precision, intent(in) :: maxAllowCourantNo ! = cfl
    double precision, intent(in) :: max_C_dx          ! = maxCelDx
      
  ! Local variables
    integer :: a
    integer :: b
    
  !-----------------------------------------------------------------------------    
  ! Calculate maximum timestep duration for numerical stability
  
    dtini = maxAllowCourantNo / max_C_dx
    a     = floor( (time - initialTime * 60.) / &
                 ( saveInterval / 60.))           
    b     = floor(((time - initialTime * 60.) + dtini / 60.) / &
                 ( saveInterval / 60.))           
    if (b > a) then
      dtini = (a + 1) * (saveInterval) - (time - initialTime * 60.) * 60.
    end if

    ! if dtini extends beyond final time, then truncate it
    if (time + dtini / 60. > tfin * 60.) dtini = (tfin * 60. - time) * 60.

  end subroutine  

  subroutine calc_dimensionless_numbers(j)
  
    IMPLICIT NONE

  !-----------------------------------------------------------------------------
  ! Description:
  !   Compute dimensionless parameters for deciding which depth computation
  !   schemes to use - normal depth or diffusive depth.
  !
  ! Method:
  !   
  !
  ! Current Code Owner: NOAA-OWP, Inland Hydraulics Team
  !
  ! Code Description:
  !   Language: Fortran 90.
  !   This code is written to JULES coding standards v1.
  !-----------------------------------------------------------------------------

  ! Subrouting arguments
    integer, intent(in) :: j
    
  ! Local variables
    integer :: i
    integer :: ncomp
    double precision :: wl_us
    double precision :: depth_us
    double precision :: q_us
    double precision :: v_us
    double precision :: pere_us
    double precision :: r_us
    double precision :: sk_us
    double precision :: ch_us
    double precision :: wl_ds
    double precision :: depth_ds
    double precision :: q_ds
    double precision :: v_ds
    double precision :: pere_ds
    double precision :: r_ds
    double precision :: sk_ds
    double precision :: ch_ds
    double precision :: ch_star_avg
    double precision :: channel_length
    double precision :: avg_celerity
    double precision :: avg_velocity
    double precision :: avg_depth
    double precision :: maxValue
    double precision :: dimlessWaveLength
      
  !-----------------------------------------------------------------------------  
    maxValue          = 1e7
    dimlessWaveLength = 4000.
    
    ncomp = frnw_g(j, 1)
    do i = 1, ncomp - 1

      ! upstream metrics
      wl_us    = newY(i, j)                          ! water level
      depth_us = newArea(i, j) / bo(i, j)            ! depth (rectangular?)
      q_us     = newQ(i, j)                          ! flow
      v_us     = abs(newQ(i, j) / newArea(i, j) )    ! velocity
      pere_us  = pere(i, j)                          ! wetted perimeter
      r_us     = newArea(i, j) / pere(i, j)          ! hydraulic radius
      sk_us    = sk(i, j)                            ! 1/Mannings' N 
      ch_us    = sk(i, j) * r_us ** (1. / 6.)        ! parameter for Chezy's constant

      ! downstream metrics
      wl_ds    = newY(i+1, j)                        ! water level
      depth_ds = newArea(i+1, j) / bo(i+1, j)        ! depth (rectangular?)
      q_ds     = newQ(i+1, j)                        ! flow
      v_ds     = abs(newQ(i+1, j) / newArea(i+1, j)) ! velocity
      pere_ds  = pere(i+1, j)                        ! wetted perimeter
      r_ds     = newArea(i+1, j) / pere(i+1, j)      ! hydraulic radius
      sk_ds    = sk(i+1, j)                          ! 1/Mannings' N 
      ch_ds    = sk(i+1, j) * r_ds ** (1./6.)        ! parameter for Chezy's constant


      ch_star_avg  = ((ch_us + ch_ds) / 2.)  / sqrt(grav)
      avg_celerity = (celerity(i, j) + celerity(i + 1, j)) / 2.0
      avg_velocity = (v_us + v_ds) / 2.
      avg_depth    = (depth_us + depth_ds) / 2.
      
      channel_length = dx(i, j)

      ! dimensionless Courant number
      dimensionless_Cr(i, j) = abs(avg_velocity / avg_celerity)
      if (dimensionless_Cr(i, j) > maxValue) dimensionless_Cr(i, j) = maxValue

      ! dimensionless Froude number
      dimensionless_Fo(i, j) = avg_velocity / sqrt(grav * avg_depth)
      if (dimensionless_Fo(i, j) > maxValue) dimensionless_Fo(i, j) = maxValue

      ! dimensionless Friction parameter (influence of friction effects on river flow)
      dimensionless_Fi(i, j) = 2. * dimensionless_Cr(i,j) / (ch_star_avg ** 2.) * &
                               dimlessWaveLength
      if (dimensionless_Fi(i, j) > maxValue) dimensionless_Fi(i, j) = maxValue

      ! dimensionless Friction parameter (influence of Fi and Courant number)
      dimensionless_Fc(i, j) = dimensionless_Cr(i, j) * dimensionless_Fi(i, j)
      if (dimensionless_Fc(i, j) > maxValue) dimensionless_Fc(i, j) = maxValue

      ! dimensionless ratio of Courant to Froude numbers or surface wave to measured wave celerity
      dimensionless_Di(i, j) = (dimensionless_Cr(i, j) / &
                                dimensionless_Fo(i, j)) ** 2.
      if (dimensionless_Di(i,j) > maxValue) dimensionless_Di(i,j) = maxValue

      ! dimensionless diffusion coefficient (ratio of wave diffusion to wave advection)
      dimensionless_D(i,j)  = dimensionless_Di(i,j) / dimensionless_Fc(i,j)
      if (dimensionless_D(i,j) .gt. maxValue) dimensionless_D(i,j) = maxValue
        
    end do
      
  end subroutine calc_dimensionless_numbers
  

  subroutine mesh_diffusive_forward(dtini_given, t0, t, tfin, saveInterval,j)

    IMPLICIT NONE
    
  !-----------------------------------------------------------------------------
  ! Description:
  !   Compute discharge using diffusive equation that is numerically solved by
  !   Crank-Nicolson + Hermite Interpolation method.
  !
  ! Method:
  !   
  ! Current Code Owner: NOAA-OWP, Inland Hydraulics Team
  !
  ! Code Description:
  !   Language: Fortran 90.
  !   This code is written to JULES coding standards v1.
  !-----------------------------------------------------------------------------

    ! Subroutine Arguments
      integer, intent(in) :: j
      double precision, intent(in) :: dtini_given
      double precision, intent(in) :: t0
      double precision, intent(in) :: t
      double precision, intent(in) :: tfin
      double precision, intent(in) :: saveInterval
    
    ! Local variables
      integer :: ncomp
      integer :: i, irow, flag_da, n
      double precision :: a1, a2, a3, a4
      double precision :: b1, b2, b3, b4
      double precision :: dd1, dd2, dd3, dd4
      double precision :: h1, h2, h3, h4
      double precision :: allqlat
      double precision :: qy, qxy, qxxy, qxxxy
      double precision :: ppi, qqi, rri, ssi, sxi
      double precision :: cour, cour2
      double precision :: alpha
      double precision :: currentQ
      double precision :: eei_ghost, ffi_ghost, exi_ghost
      double precision :: fxi_ghost, qp_ghost, qpx_ghost
    !-----------------------------------------------------------------------------
    !* change 20210228: All qlat to a river reach is applied to the u/s boundary
    !* Note: lateralFlow(1,j) is already added to the boundary
    
      eei = -999.
      ffi = -999.
      exi = -999.
      fxi = -999.
      eei(1) = 1.
      ffi(1) = 0.
      exi(1) = 0.
      fxi(1) = 0.

      ncomp  = frnw_g(j, 1)
    
      ! sum of lateral inflows along the reach
      !allqlat = sum(lateralFlow(2:ncomp - 1, j) * dx(2:ncomp - 1, j))
      allqlat = 0.0
      do i = 2, ncomp - 1
        allqlat = allqlat + lateralFlow(i, j) * dx(i, j)
      end do

      !print *, '****** DIFFUSIVE FORWARD *****'
      !print *, '---------'
      ncomp = frnw_g(j, 1)
      
      do i = 2, ncomp          
        
        cour  = dtini / dx(i - 1, j)
        cour2 = abs(celerity(i, j) ) * cour
               
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

        qy    = a1 * oldQ(i-1,j) + a2 * oldQ(i,j) + &
                a3 * qpx(i-1,j) + a4 * qpx(i,j)
        qxy   = b1 * oldQ(i-1,j) + b2 * oldQ(i,j) + &
                b3 * qpx(i-1,j) + b4 * qpx(i,j)
        qxxy  = dd1* oldQ(i-1,j) + dd2* oldQ(i,j) + &
                dd3* qpx(i-1,j) + dd4* qpx(i,j)
        qxxxy = h1 * oldQ(i-1,j) + h2 * oldQ(i,j) + &
                h3 * qpx(i-1,j) + h4 * qpx(i,j)

        ppi = - theta * diffusivity(i,j) * dtini / &
                ( dx(i-1,j) ** 2.0 ) * 2.0 / (alpha*(alpha + 1.0)) * alpha
        qqi = 1.0 - ppi * (alpha + 1.0) / alpha
        rri = ppi / alpha

        ssi = qy  + dtini * diffusivity(i,j) * ( 1.0 - theta ) * qxxy
        sxi = qxy + dtini * diffusivity(i,j) * ( 1.0 - theta ) * qxxxy

        eei(i) = -1.0 * rri / ( ppi * eei(i-1) + qqi )
        ffi(i) = ( ssi - ppi * ffi(i-1) ) / ( ppi * eei(i-1) + qqi )
        exi(i) = -1.0 * rri / ( ppi * exi(i-1) + qqi )
        fxi(i) = ( sxi - ppi * fxi(i-1) ) / ( ppi * exi(i-1) + qqi )
        
      end do
      
      ! Ghost point calculation
      cour  = dtini / dx(ncomp-1, j)
      cour2 = abs(celerity(ncomp-1, j)) * cour

      a1 = 3.0 * cour2 ** 2.0 - 2.0 * cour2 ** 3.0
      a2 = 1 - a1
      a3 = ( cour2 ** 2.0 - cour2 ** 3.0 ) * dx(ncomp-1, j)
      a4 = ( -1.0 * cour2 + 2.0 * cour2 ** 2.0 - cour2 ** 3.0 ) * dx(ncomp-1, j)

      b1 = ( 6.0 * cour2 - 6.0 * cour2 ** 2.0 ) / ( -1.0 * dx(ncomp-1, j) )
      b2 = - b1
      b3 = ( 2.0 * cour2 - 3.0 * cour2 ** 2.0 ) * ( -1.0 )
      b4 = ( -1.0 + 4.0 * cour2 - 3.0 * cour2 ** 2.0 ) * ( -1.0 )

      dd1 = ( 6.0 - 12.0 * cour2 ) / ( dx(ncomp-1, j) ** 2.0 )
      dd2 = - dd1
      dd3 = ( 2.0 - 6.0 * cour2 ) / dx(ncomp-1, j)
      dd4 = ( 4.0 - 6.0 * cour2 ) / dx(ncomp-1, j)

      h1 = 12.0 / ( dx(ncomp-1, j) ** 3.0 )
      h2 = - h1
      h3 = 6.0 / ( dx(ncomp-1, j) ** 2.0 )
      h4 = h3

      alpha = 1.0

      qy    = a1 * oldQ(ncomp, j) + a2 * oldQ(ncomp-1, j) + &
              a3 * qpx(ncomp, j) + a4 * qpx(ncomp-1, j)
      qxy   = b1 * oldQ(ncomp, j) + b2 * oldQ(ncomp-1, j) + &
              b3 * qpx(ncomp, j) + b4 * qpx(ncomp-1, j)
      qxxy  = dd1 * oldQ(ncomp, j) + dd2 * oldQ(ncomp-1, j) + &
              dd3 * qpx(ncomp, j) + dd4 * qpx(ncomp-1, j)
      qxxxy = h1 * oldQ(ncomp, j) + h2 * oldQ(ncomp-1, j) + &
              h3 * qpx(ncomp, j) + h4 * qpx(ncomp-1, j)

      ppi = - theta * diffusivity(ncomp, j) * dtini / &
              ( dx(ncomp-1, j) ** 2.0 ) * 2.0 / (alpha*(alpha + 1.0)) * alpha
      qqi = 1.0 - ppi * (alpha + 1.0) / alpha
      rri = ppi / alpha

      ssi = qy  + dtini * diffusivity(ncomp-1, j) * ( 1.0 - theta ) * qxxy
      sxi = qxy + dtini * diffusivity(ncomp-1, j) * ( 1.0 - theta ) * qxxxy

      eei_ghost = -1.0 * rri / ( ppi * eei(ncomp) + qqi )
      ffi_ghost = ( ssi - ppi * ffi(ncomp) ) / ( ppi * eei(ncomp) + qqi )

      exi_ghost = -1.0 * rri / ( ppi * exi(ncomp) + qqi )
      fxi_ghost = ( sxi - ppi * fxi(ncomp) ) / ( ppi * exi(ncomp) + qqi )

      qp_ghost  = oldQ(ncomp-1, j)
      qpx_ghost = 0.
      
      ! when a reach has usgs streamflow data at its location, apply DA
      !if (usgs_da_reach(j) /= 0) then
        
      !  allocate(varr_da(nts_da))
      !  do n = 1, nts_da
      !      varr_da(n) = usgs_da(n, j)
      !  end do
      !  qp(ncomp,j) = intp_y(nts_da, tarr_da, varr_da, t + dtini/60.)
      !  flag_da = 1
        ! check usgs_da value is in good quality
      !  irow = locate(tarr_da, t + dtini/60.)
      !  if (irow == nts_da) then
      !    irow = irow-1
      !  endif        
      !  if ((varr_da(irow)<= -4443.999).or.(varr_da(irow+1)<= -4443.999)) then
          ! when usgs data is missing or in poor quality
          qp(ncomp,j)  = eei(ncomp) * qp_ghost + ffi(ncomp)
          flag_da = 0
      !  endif
      !  deallocate(varr_da)
     ! else
        
        qp(ncomp,j)  = eei(ncomp) * qp_ghost + ffi(ncomp)
        flag_da = 0
     ! endif
      
      qpx(ncomp,j) = exi(ncomp) *qpx_ghost + fxi(ncomp)

      do i = ncomp-1, 1, -1
        qp(i, j)  = eei(i) * qp(i+1, j) + ffi(i)
        qpx(i, j) = exi(i) * qpx(i+1, j) + fxi(i)
     end do

      ! when a reach hasn't been applied to DA 
     ! if ((usgs_da_reach(j) == 0).or.(flag_da == 0)) then
       qp(1, j) = newQ(1, j)
       qp(1, j) = qp(1, j) + allqlat
     ! endif
    
      do i = 1, ncomp
        if (abs(qp(i, j)) < q_llm) then
          qp(i, j) = q_llm    
        end if
      end do

      ! update newQ
      do i = 1, ncomp
        newQ(i, j) = qp(i, j)
      end do
    ! ============== DEBUG to find unstable flow calcls =================
    !        do i= ncomp, 1, -1
    !            if (abs(qp(i,j)) .gt. 2E4) then
    !                print *, 'j:', j
    !                print *, 'i:', i
    !                print *, 'ncomp', ncomp
    !                print *, 't:', t
    !                print *, 'calculated flow value:', qp(i,j)
    !                print *, 'previous upstream flow', oldQ(i-1,j)
    !                print *, 'diffusivity(i,j):', diffusivity(i,j)
    !                print *, 'eei(i):', eei(i)
    !                print *, 'qp(i+1,j):', qp(i+1,j)
    !                print *, 'ffi(i):', ffi(i)
    !                print *, 'qp_ghost:', qp_ghost
    !                print *, 'ffi(ncomp):', ffi(ncomp)
    !                print *, 'eei(ncomp):', eei(ncomp)
    !                print *, 'cour2:', cour2
    !                print *, 'qp(ncomp,j):', qp(ncomp,j)
    !                print *, 'allqlat:', allqlat
    !                print *, 'upstream inflow newQ(1,j):', newQ(1,j)
    !                stop
    !            end if
    !        end do
      
  end subroutine mesh_diffusive_forward
  
  subroutine mesh_diffusive_backward(dtini_given, t0, t, tfin, saveInterval, j) ! rightBank)

    IMPLICIT NONE

  !-----------------------------------------------------------------------------
  ! Description:
  !   Compute water depth using diffusive momentum equation or normal depth
  !   with computed Q from the forward step.
  !
  ! Method:
  !   
  !
  ! Current Code Owner: NOAA-OWP, Inland Hydraulics Team
  !
  ! Code Description:
  !   Language: Fortran 90.
  !   This code is written to JULES coding standards v1.
  !-----------------------------------------------------------------------------

  ! Subroutine arguments
    integer, intent(in) :: j
    double precision, intent(in) :: dtini_given
    double precision, intent(in) :: t0
    double precision, intent(in) :: t
    double precision, intent(in) :: tfin
    double precision, intent(in) :: saveInterval
      
  ! Subroutine variables
    integer :: depthCalOk(mxncomp)
    integer :: i, ncomp  
    integer :: jj
    double precision :: xt
    double precision :: q_sk_multi, sfi
    double precision :: vel, currentQ
    double precision :: S_ncomp                             
    double precision :: tempDepthi_1
    double precision :: Q_cur, Q_ds, z_cur, z_ds, y_cur, y_ds
    double precision :: C_ulm

  !-----------------------------------------------------------------------------
    ncomp = frnw_g(j, 1)
    S_ncomp = (-z(ncomp, j) + z(ncomp-1, j)) / dx(ncomp-1, j)
    depthCalOk(ncomp) = 1
    elevTable = xsec_tab(1, :,ncomp, j)
    areaTable = xsec_tab(2, :,ncomp, j)
    topwTable = xsec_tab(6, :,ncomp, j)
    
    ! Estimate channel area @ downstream boundary
    call r_interpol(elevTable, areaTable, nel, &
                    newY(ncomp, j), newArea(ncomp, j))
    if (newArea(ncomp,j) .eq. -9999) then
      print*, 'At j = ',j,'i = ',ncomp, 'time =',t, 'newY(ncomp, j)=', newY(ncomp, j), &
            'newArea(ncomp, j)=', newArea(ncomp,j), 'interpolation of newArea was not possible'
!     stop
    end if

    ! Estimate channel bottom width @ downstream boundary
    call r_interpol(elevTable, topwTable, nel, &
                    newY(ncomp, j), bo(ncomp, j))
    if (bo(ncomp, j) .eq. -9999) then
      print*, 'At j = ',j,', i = ',ncomp, 'time =',t, &
              'interpolation of bo was not possible'
!     stop
    end if

    do i = ncomp, 1, -1 ! loop through reach nodes [downstream-upstream]
      currentQ = qp(i, j)
      q_sk_multi=1.0

      elevTable = xsec_tab(1,:,i,j)
      convTable = xsec_tab(5,:,i,j)
      areaTable = xsec_tab(2,:,i,j)
      pereTable = xsec_tab(3,:,i,j)
      topwTable = xsec_tab(6,:,i,j)
      skkkTable = xsec_tab(11,:,i,j)
              
      xt=newY(i, j)      
 
      ! Estimate co(i) (???? what is this?) by interpolation
      currentSquareDepth = (elevTable - z(i, j)) ** 2.
      call r_interpol(currentSquareDepth, convTable, nel, &
                      (newY(i, j)-z(i, j)) ** 2.0, co(i)) 
                   
      if (co(i) .eq. -9999) then
        ! test                
        print*, t, i, j, newY(i,j), newY(i, j)-z(i, j), co(i)  
        print*, 'At j = ',j,', i = ',i, 'time =',t, &
                'interpolation of conveyance was not possible, wl', &
                newY(i,j), 'z',z(i,j),'previous wl',newY(i+1,j) !, &
                !'previous z',z(i+1,j), 'dimensionless_D(i,j)', &
                !dimensionless_D(i,j)
  !      stop
         pause !test
      end if
      co(i) =q_sk_multi * co(i)

      ! Estimate channel area by interpolation
      call r_interpol(elevTable, areaTable, nel, &
                      xt, newArea(i,j))
      if (newArea(i,j) .eq. -9999) then
        print*, 'At j = ',j,', i = ',i, 'time =',t, &
                'interpolation of newArea was not possible'
!        stop
      end if
      
      ! Estimate wetted perimeter by interpolation
      call r_interpol(elevTable, pereTable, nel, &
                      xt, pere(i, j))
      
      ! Estimate width by interpolation
      call r_interpol(elevTable, topwTable, nel, &
                      xt, bo(i, j))
      
      ! Estimate sk(i,j) by interpolation (???? what is sk?)
      call r_interpol(elevTable, skkkTable, nel, &
                      xt, sk(i, j))

      sfi = qp(i, j) * abs(qp(i, j)) / (co(i) ** 2.0)

      ! ???? What exactly is happening, here ????
      !if (depthCalOk(i) .eq. 1) then
      
        ! Calculate celerity, diffusivity and velocity
        celerity2(i)    = 5.0 / 3.0 * abs(sfi) ** 0.3 * &
                          abs(qp(i, j)) ** 0.4 / bo(i, j) ** 0.4 &
                          / (1. / (sk(i, j) * q_sk_multi)) ** 0.6

        ! put upper limit on computed celerity to make sure internal temporal increment not being too small.
        if (i.gt.1) then
          C_ulm = cfl*dx(i-1,j)/dtini_min
        else 
          C_ulm = cfl*dx(i,j)/dtini_min
        endif

        if (celerity2(i).gt.C_ulm) then
          celerity2(i) = C_ulm
        endif
        
        diffusivity2(i) = abs(qp(i, j)) / 2.0 / bo(i, j) / abs(sfi)
        vel             = qp(i, j) / newArea(i, j)
        
        ! Check celerity value
        !if (celerity2(i) .gt. 3.0 * vel) celerity2(i) = vel * 3.0
      !else
      !  if (qp(i, j) .lt. 1) then
      !    celerity2(i) = C_llm
      !  else
      !    celerity2(i) = 1.0
      !  end if
      !  diffusivity2(i) =  diffusivity(i,j)
      !end if
      
      if (i .gt. 1) then       
        ! ====================================================
        ! Combination of Newton-Raphson and Bisection
        ! ====================================================                    
        Q_cur        = qp(i-1, j)
        Q_ds         = qp(i, j)
        z_cur        = z(i-1, j)
        z_ds         = z(i, j)
        y_ds         = newY(i, j) - z(i, j)
        y_ds         = max(y_ds, 0.005)       
        y_cur        = rtsafe(i-1, j, Q_cur, Q_ds, z_cur, z_ds, y_ds)
        tempDepthi_1 = y_cur              
        newY(i-1,j)  = tempDepthi_1 + z(i-1, j)

        if (newY(i-1, j) > 10.0**5.) newY(i-1, j) = 10.0**5.
              
        if (newY(i - 1, j) - z(i-1, j) <= 0.) then
            print *, ' newY(i-1,j)-z(i-1,j): ', newY(i-1,j)-z(i-1,j)
            print *, ' newY(i-1,j): ', newY(i-1,j)
            print *, 'z(i-1,j): ', z(i-1,j)
            print *, 'qp(i-1,j): ', qp(i-1,j)
            print*, 'depth is negative at time=,', t,'j= ', j,'i=',i-1, &
                    'newY =', (newY(jj,j),jj=1,ncomp)
            !print*, 'dimensionless_D',(dimensionless_D(jj,j),jj=1,ncomp)
            print*, 'newQ',(newQ(jj,j),jj=1,ncomp)
            print*, 'Bed',(z(jj,j),jj=1,ncomp)
            print*, 'dx',(dx(jj,j),jj=1,ncomp-1)
        end if     
      end if
           
    end do

    celerity(1:ncomp, j) =  sum(celerity2(1:ncomp)) / ncomp
    
    if (celerity(1, j) < C_llm) celerity(1:ncomp,j) = C_llm
    
    diffusivity(1:ncomp, j) = sum(diffusivity2(1:ncomp)) / ncomp
    
    do i = 1, ncomp
      if (diffusivity(i, j) > D_ulm) diffusivity(i, j) = D_ulm !!! Test
      if (diffusivity(i, j) < D_llm) diffusivity(i, j) = D_llm !!! Test
    end do
  end subroutine mesh_diffusive_backward

  function rtsafe(i, j, Q_cur, Q_ds, z_cur, z_ds, y_ds)
    
    implicit none
        
  !----------------------------------------------------------------------------------------------------------
  ! Description:
  !   Compute water depth using diffusive momentum equation using a combination of Newton-Raphson and Bisection
  !
  ! Method:     
  !  p.1189, Numerical Recipes in F90  
  !  Using a combination of newton-raphson and bisection, find the root of a function bracketed
  !  between x1 and x2. the root, returned as the function value rtsafe, will be refined until
  !  its accuracy is known within ±xacc.
  !  - funcd is a user-supplied subroutine that returns both the function value and the first   
  !    derivative of the function.
  !  - parameter: maxit is the maximum allowed number of iterations.
  !----------------------------------------------------------------------------------------------------------

    ! Subroutine arguments
    integer,          intent(in) :: i, j
    double precision, intent(in) :: Q_cur, Q_ds, z_cur, z_ds, y_ds
  
    ! Subroutine local variable
    integer,          parameter :: maxit = 40
    double precision, parameter :: xacc = 1e-4
    integer                     :: iter
    integer                     :: xcolID, ycolID
    double precision            :: x1, x2, df, dxx, dxold, f, fh, fl, temp, xh, xl
    double precision            :: y_norm, y_ulm_multi, y_llm_multi, elv_norm, y_old
    double precision            :: rtsafe    
    
    y_ulm_multi = 2.0
    y_llm_multi = 0.1
    
    xcolID   = 10
    ycolID   = 1
    elv_norm = intp_xsec_tab(i, j, nel, xcolID, ycolID, abs(Q_cur)) ! normal elevation not depth
    y_norm   = elv_norm - z(i, j)
    y_old    = oldY(i, j) - z(i, j)    
   
    ! option 1 for initial point
    !x1       = y_norm * y_llm_multi
    !x2       = y_norm * y_ulm_multi
    ! option 2 for initial point
    x1       = 0.5 * (y_norm + y_old) * y_llm_multi
    x2       = 0.5 * (y_norm + y_old) * y_ulm_multi
       
    call funcd_diffdepth(i, j, Q_cur, Q_ds, z_cur, z_ds, x1, y_ds, fl, df)

    call funcd_diffdepth(i, j, Q_cur, Q_ds, z_cur, z_ds, x2, y_ds, fh, df)    
   
    if ((fl > 0.0 .and. fh > 0.0) .or. (fl < 0.0 .and. fh < 0.0)) then
      rtsafe = y_norm
      return
    endif

    if (fl == 0.0) then
      rtsafe = x1
      return
    elseif (fh == 0.0) then
      rtsafe = x2
      return
    elseif (fl < 0.0) then ! orient the search so that f(xl) < 0.
      xl = x1
      xh = x2
    else
      xh = x1
      xl = x2
    end if

    rtsafe = 0.50 * (x1 + x2)      ! initialize the guess for root
    dxold  = abs(x2 - x1)          ! the “stepsize before last,”
    dxx    = dxold                 ! and the last step.
    
    call funcd_diffdepth(i, j, Q_cur, Q_ds, z_cur, z_ds, rtsafe, y_ds, f, df)

    do iter = 1, maxit             ! loop over allowed iterations.
      if (((rtsafe - xh) * df - f) * ((rtsafe - xl) * df - f) > 0.0 .or. &
                                    abs(2.0 * f) > abs(dxold * df) ) then
      ! bisect if newton out of range, or not decreasing fast enough.
        dxold  = dxx
        dxx    = 0.50 * (xh - xl)
        rtsafe = xl + dxx
        if (xl == rtsafe) return       
      else                         ! newton step acceptable. take it.
        dxold  = dxx
        dxx    = f / df
        temp   = rtsafe
        rtsafe = rtsafe - dxx
        if (temp == rtsafe) return      
      end if
      
      if (abs(dxx) < xacc) return 

      ! one new function evaluation per iteration.
      call funcd_diffdepth(i, j, Q_cur, Q_ds, z_cur, z_ds, rtsafe, y_ds, f, df)

      if (f < 0.0) then            ! maintain the bracket on the root.
        xl = rtsafe
      else
        xh = rtsafe
      end if
    end do

    ! when root is not converged:
    rtsafe = y_norm
    
  end function rtsafe

  subroutine funcd_diffdepth(i, j, Q_cur, Q_ds, z_cur, z_ds, y_cur, y_ds, f, df)
    
    implicit none
    
    !-------------------------------------------------------------------------
    ! Description:
    !   Compute diffusive momentum function value and the first derivative
    !
    ! Method:     
    !   Analytical function and its analytical derivative
    !-------------------------------------------------------------------------

    ! subroutine arguments
    integer,          intent(in)  :: i, j
    double precision, intent(in)  :: Q_cur, Q_ds, z_cur, z_ds, y_cur, y_ds
    double precision, intent(out) :: f, df
    
    ! subroutine local variables
    integer          :: xcolID, ycolID
    double precision :: elv_cur, elv_ds, conv_cur, conv_ds, sf_cur, sf_ds, slope
    double precision :: dKdA_cur, topw_cur

    xcolID  = 1
    ! f(y_cur): function of y at the current node
    ! - energy slope at downstream node
    ycolID  = 5
    elv_ds  = y_ds + z_ds
    conv_ds = intp_xsec_tab(i + 1, j, nel, xcolID, ycolID, elv_ds)
    sf_ds   = abs(Q_ds) * Q_ds / conv_ds**2.0    
    ! - energy slope at current node
    elv_cur  = y_cur + z_cur
    conv_cur = intp_xsec_tab(i, j, nel, xcolID, ycolID, elv_cur)
    sf_cur   = abs(Q_cur) * Q_cur / conv_cur**2.0
    ! - f(y_cur)
    slope    = (z(i, j) - z(i+1, j)) / dx(i, j)
    slope    = max(slope, so_llm)
    f        = y_cur - y_ds + slope * dx(i, j) - 0.50 * (sf_cur + sf_ds) * dx(i, j)

    ! df/dy at y at current node
    ! - dK/dA
    ycolID   = 9
    dKdA_cur = intp_xsec_tab(i, j, nel, xcolID, ycolID, elv_cur)
    ! - top width
    ycolID   = 6
    topw_cur = intp_xsec_tab(i, j, nel, xcolID, ycolID, elv_cur)
    df = 1.0 + (abs(Q_cur) * Q_cur / conv_cur**3.0) * dx(i, j) * topw_cur * dKdA_cur

  end subroutine funcd_diffdepth

  double precision function intp_xsec_tab(i, j, nrow, xcolID, ycolID, x)
    
    implicit none
    
    !-------------------------------------------------------------------
    ! Description:
    !   Interpolate given columns of hydraulic lookup table
    !
    ! Method:     
    !   linear interpolation between selected adjacent data points
    !-------------------------------------------------------------------

    ! subroutine arguments
    integer         , intent(in) :: i, j, nrow, xcolID, ycolID
    double precision, intent(in) :: x
    
    ! subroutine local variables
    integer                           :: irow
    double precision                  :: x1, y1, x2, y2, y
    double precision, dimension(nrow) :: xarr, yarr

    xarr = xsec_tab(xcolID, 1:nrow, i, j)
    yarr = xsec_tab(ycolID, 1:nrow, i, j)

    irow = locate(xarr, x)

    if (irow == 0)    irow = 1
    if (irow == nrow) irow = nrow-1

    x1            =   xarr(irow)
    y1            =   yarr(irow)
    x2            =   xarr(irow+1)
    y2            =   yarr(irow+1)
    y             = LInterpol(x1, y1, x2, y2, x)
    intp_xsec_tab = y
  end function intp_xsec_tab

  subroutine r_interpol(x, y, kk, xrt, yt)
      
      implicit none
      !---------------------------------------------------------------------------
      ! Description:
      !   Estimate y for a given along x and y arrays using linear interpolation
      !---------------------------------------------------------------------------         
      ! subroutine arguments
      integer,          intent(in)  :: kk
      double precision, intent(in)  :: xrt, x(kk), y(kk)
      double precision, intent(out) :: yt
      ! subroutine local variables
      integer :: k

      if ((xrt <= maxval(x)) .and. (xrt >= minval(x))) then
      
        do k = 1, kk-1
          if (((x(k) - xrt) * (x(k+1) - xrt)) <= 0.0) then
            yt = (xrt - x(k)) / (x(k+1) - x(k)) * (y(k+1) - y(k)) + y(k)
            EXIT
          endif
        end do
      
      else if (xrt >= maxval(x)) then
!       print*, xrt, ' the given x data point is larger than the upper limit of the set of x data points'
!       print*, 'the upper limit: ', maxval(x)
        yt = (xrt - x(kk-1)) / (x(kk) - x(kk-1)) * (y(kk) - y(kk-1)) + y(kk-1) ! extrapolation

      else
!       print*, xrt, ' the given x data point is less than lower limit of the range of known x data point, '
!       print*, 'so linear interpolation cannot be performed.'
!        yt = -9999.0
        yt = minval(y)
!       print*, 'The proper range of x is that: ', 'the upper limit: ', maxval(x),&
!               ' and lower limit: ', minval(x)
!       print*, 'kk', kk
!       print*, 't', 'i', dmyi, 'j', dmyj
!       print*, 'x', (x(k), k=1, kk)
!       print*, 'y', (y(k), k=1, kk)
      end if
    
  end subroutine r_interpol

  subroutine normal_crit_y(i, j, q_sk_multi, So, dsc, y_norm, y_crit, area_n, area_c)
      
    implicit none
        
    !-------------------------------------------------------------------------------------------------
    ! Description:
    !   Estimate normal and critical depth and related areas
    !
    ! Method:
    !   normal depth by linearly interpolating elevation & conveyance columns w.r.t. given conveyance
    !   normal depth area by linearly interpolation elevation & area w.r.t. given normal depth
    !   critical depth by iterative method and the computed depth leads to the estimating of the area	
    !--------------------------------------------------------------------------------------------------

    ! subroutine arguments
    integer,          intent(in)  :: i, j
    double precision, intent(in)  :: q_sk_multi, So, dsc
    double precision, intent(out) :: y_norm, y_crit, area_n, area_c
      
    ! subroutine local variables 
    double precision :: area_0, width_0, errorY

    elevTable = xsec_tab(1, :, i, j)
    areaTable = xsec_tab(2, :, i, j)
    pereTable = xsec_tab(3, :, i, j)
    convTable = xsec_tab(5, :, i, j)
    topwTable = xsec_tab(6, :, i, j)
      
    call r_interpol(convTable, areaTable, nel, dsc / sqrt(So), area_n)
    call r_interpol(convTable, elevTable, nel, dsc / sqrt(So), y_norm)
    call r_interpol(elevTable, areaTable, nel, oldY(i, j), area_0)  ! initial estimate for critical depth
    call r_interpol(elevTable, topwTable, nel, oldY(i, j), width_0) ! initial estimate for critical depth

    area_c = area_0
    errorY = 100.
       
    do while (errorY > 0.0001)
      area_c = (dsc * dsc * width_0 / grav) ** (1./3.)
      errorY = abs(area_c - area_0)
        
      call r_interpol(areaTable, topwTable, nel, area_c, width_0)
        
      area_0 = area_c
    end do
    
    call r_interpol(areaTable,elevTable,nel,area_c,y_crit)
      
    if (y_norm .eq. -9999) then
      print*, 'At j = ',j,', i = ',i, 'interpolation of y_norm in calculating normal area was not possible, Q', &
              dsc,'slope',So
!     stop
    end if
  end subroutine normal_crit_y

  double precision function LInterpol(x1, y1, x2, y2, x)
      
    implicit none
      
    !-------------------------------------------------------------------------------------
    ! Description:
    !   Estimate y for a given x applying linear interpolation to two given (x, y) points
    !-------------------------------------------------------------------------------------          
      
    ! function arguments
    double precision, intent(in) :: x1, y1, x2, y2, x
      
    if (abs(x2-x1).lt.0.0001) then
      ! to prevent absurdly small value in the denominator
      LInterpol = 0.5*(y1 + y2)
    else
      LInterpol = (y2 - y1) / (x2 - x1) * (x - x1) + y1
    endif
      
  end function LInterpol

  double precision function intp_y(nrow, xarr, yarr, x)
      
    implicit none
    !-------------------------------------------------------------------------------------
    ! Description:
    !   Estimate y for a given x applying linear interpolation to two x and y arrays
    !-------------------------------------------------------------------------------------         

    ! function arguments
    integer,                           intent(in) :: nrow
    double precision,                  intent(in) :: x      
    double precision, dimension(nrow), intent(in) :: xarr, yarr
    ! function local variables
    integer          :: irow
    double precision :: x1, y1, x2, y2, y

    irow = locate(xarr, x)
      
    if (irow == 0)    irow = 1
    if (irow == nrow) irow = nrow - 1
      
    x1     = xarr(irow)
    y1     = yarr(irow)
    x2     = xarr(irow+1)
    y2     = yarr(irow+1)
    y      = LInterpol(x1, y1, x2, y2, x)
    intp_y = y

  end function intp_y

  integer function locate(xx, x)
      
    implicit none
      
    !------------------------------------------------------------------------------------
    ! Description:               
    !   Locate function in p.1045,Numerical Recipes in Fortran f90
    !      
    ! Method:
    !   klo=max(min(locate(xa,x),n-1),1) In the Fortran 77 version of splint,
    !   there is in-line code to find the location in the table by bisection. Here
    !   we prefer an explicit call to locate, which performs the bisection. On
    !   some massively multiprocessor (MMP) machines, one might substitute a different,
    !   more parallel algorithm (see next note).
    !  Given an array xx(1:N), and given a value x, returns a value j such that x is between
    !  xx(j) and xx(j + 1). xx must be monotonic, either increasing or decreasing.
    !  j = 0 or j = N is returned to indicate that x is out of range.
    !------------------------------------------------------------------------------------
      
    ! function arguments
    double precision,               intent(in) :: x      
    double precision, dimension(:), intent(in) :: xx

    ! function local variables
    integer :: n, jl, jm, ju
    logical :: ascnd

    n     = size(xx)
    ascnd = (xx(n) >= xx(1))    ! True if ascending order of table, false otherwise.
    jl    = 0                   ! Initialize lower
    ju    = n + 1               ! and upper limits.
      
    do
      if ((ju - jl) <= 1) exit  ! Repeat until this condition is satisfied.
        
      jm = (ju + jl) / 2        ! Compute a midpoint,
        
      if (ascnd .eqv. (x >= xx(jm))) then
        jl = jm                 ! and replace either the lower limit
      else
        ju = jm                 ! or the upper limit, as appropriate.
      end if
    end do

    if (x == xx(1)) then        ! Then set the output, being careful with the endpoints.
      locate = 1
    else if (x == xx(n)) then
      locate = n - 1
    else
      locate = jl
    end if

  end function locate

end module diffusive
