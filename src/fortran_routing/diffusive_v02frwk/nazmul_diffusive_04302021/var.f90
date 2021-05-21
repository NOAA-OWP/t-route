module var

    implicit none
    save

    real, parameter :: grav = 9.81
    real, parameter :: TOLERANCE = 1e-8
    real, allocatable :: area(:), bo(:,:)
    real, allocatable :: depth(:)
    real, allocatable :: co(:), dbdx(:)
    real, allocatable :: dx(:,:), froud(:), courant(:)
    real, allocatable :: dt(:)
	!**arrays for branching channel application
    integer, allocatable :: ndep(:), uslinks(:,:), dslink(:), instrdflag(:,:)
    real, allocatable :: areap(:, :), qp(:, :), z(:, :), sk(:, :)
    real, allocatable :: dqp(:,:), dap(:,:), dqc(:,:), dac(:,:)
    real, allocatable :: celerity(:,:), diffusivity(:,:), diffusivity2(:), celerity2(:)
    real, allocatable :: eei(:), ffi(:), exi(:), fxi(:), qpx(:,:), qcx(:)
    real, allocatable :: USBoundary(:,:,:), DSBoundary(:,:,:)
    integer, allocatable :: upBoundTableEntry(:), downBoundTableEntry(:)
    ! change for unsteady flow
    real, allocatable :: pere(:,:),dpda(:)
    real, allocatable :: oldQ(:,:), newQ(:,:), oldArea(:,:), newArea(:,:), oldY(:,:), newY(:,:)
    real, allocatable :: lateralFlow(:,:)
    real, allocatable :: dimensionless_Cr(:,:), dimensionless_Fo(:,:), dimensionless_Fi(:,:)
    real, allocatable :: dimensionless_Di(:,:), dimensionless_Fc(:,:), dimensionless_D(:,:)
    real, allocatable :: ini_y(:), ini_q(:)
    integer, allocatable ::  noQSKtable(:)
    real, allocatable :: lowerLimitCount(:), higherLimitCount(:), volRemain(:,:)
    integer, allocatable :: currentROutingDiffusive(:), notSwitchRouting(:)
    real :: minDx, maxCelerity
    integer, allocatable :: currentRoutingNormal(:,:), routingNotChanged(:,:)
    real(kind=4), allocatable :: elevTable(:),areaTable(:), skkkTable(:)
    real(kind=4), allocatable :: pereTable(:),rediTable(:)
    real(kind=4), allocatable :: convTable(:),topwTable(:)
    real(kind=4), allocatable :: nwi1Table(:),dPdATable(:)
    real(kind=4), allocatable :: ncompElevTable(:), ncompAreaTable(:)
    real(kind=4), allocatable :: xsec_tab(:, :, :, :)
    real(kind=4), allocatable :: currentSquareDepth(:)
    integer :: maxTableLength, nel
    character*4 :: file_num
    integer :: ots, option_dsbc
    real :: cfl, qus, f, yds,  us, thes
    real :: vv
    real :: dtini
    real :: frus2, S_0, y_norm_us, y_crit_us, area_norm, area_crit, minNotSwitchRouting, minNotSwitchRouting2
    integer :: ncomp
    !* variables for branching channels
    integer :: nlinks
    real :: theta
    real :: dxini,lastKnownDiffuDT

    real, dimension(:,:), allocatable :: z_ar_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g, so_ar_g
    real, dimension(:,:), allocatable :: mann_ar_g, manncc_ar_g, dx_ar_g
    integer, dimension(:,:), allocatable :: frnw_g
    real, dimension(:,:), allocatable :: dfrnw_g
    real(KIND=4):: z_g, bo_g, traps_g, tw_g, twcc_g, so_g, mann_g, manncc_g
    real, dimension(:,:), allocatable :: ubcd_g
    real, dimension(:), allocatable ::  dbcd_g
    real, dimension(:,:,:), allocatable :: qlat_g
    integer :: applyNaturalSection
    integer :: mxncomp_g
    real :: cfl_g, dt_ub_g, dt_db_g, dt_ql_g, dtini_g
    real :: t0_g, tfin_g, tfin_qlat_g
    real :: saveInterval_g, saveInterval_ev_g
    integer :: ntss_g, nts_ql_g, nts_ub_g, nts_db_g, ntss_ev_g, nl_ubcd_ar_g
    integer :: nel_g, nrch_g
    real :: timesdepth_g
endmodule var
