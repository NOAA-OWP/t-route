module arrays_module

    implicit none
    save

    double precision, allocatable :: area(:), bo(:,:)
    double precision, allocatable :: av11(:), av12(:), av21(:), av22(:)
    double precision, allocatable ::  ci1(:), ci2(:)
    double precision, allocatable :: aso(:,:), f1(:), f2(:), depth(:)
    double precision, allocatable :: g11inv(:), g12inv(:), g21inv(:), g22inv(:)
    double precision, allocatable :: b11(:), b12(:), b21(:), b22(:)
    double precision, allocatable :: eps2(:), eps4(:), d1(:), d2(:), u(:), c(:)
    double precision, allocatable :: co(:), gso(:,:), dbdx(:,:)
    double precision, allocatable :: dx(:,:), froud(:), courant(:)
    double precision, allocatable :: dt(:)

	!**arrays for branching channel application
    integer, allocatable :: ndep(:), uslinks(:,:), dslink(:), instrdflag(:,:), nx1(:)
    double precision, allocatable :: areap(:, :), qp(:, :, :), z(:, :), sk(:, :)            ! change 20210628
    double precision, allocatable :: dqp(:,:), dap(:,:), dqc(:,:), dac(:,:)
    double precision, allocatable :: celerity(:,:),velocity(:,:), diffusivity(:,:), diffusivity2(:), celerity2(:)
    double precision, allocatable :: eei(:), ffi(:), exi(:), fxi(:), qpx(:,:), qcx(:)
    double precision, allocatable :: USBoundary(:,:,:), DSBoundary(:,:,:)
    integer, allocatable :: upBoundTableEntry(:), downBoundTableEntry(:)
    ! change for unsteady flow
    double precision, allocatable :: pere(:,:),dpda(:)
    double precision, allocatable :: oldQ(:,:), newQ(:,:,:), oldArea(:,:), newArea(:,:)
	double precision, allocatable :: oldY(:,:), newY(:,:), normalDepthAtNodes(:,:)      ! change 20210628
    double precision, allocatable :: added_Q(:,:,:)             !change 20210713

    integer, allocatable :: ityp(:), latFlowLocations(:,:), dataInEachLatFlow(:,:), latFlowType(:,:), latFlowXsecs(:,:)
    double precision, allocatable :: lateralFlowTable(:,:,:,:), lateralFlow(:,:,:)
    ! for additional lateral flow of the structures
    integer, allocatable :: latFlowLocations2(:,:), dataInEachLatFlow2(:,:), latFlowType2(:,:), latFlowXsecs2(:,:), noLatFlow2(:)
    double precision, allocatable :: lateralFlowTable2(:,:,:,:), lateralFlow2(:,:)

    double precision, allocatable :: dimensionless_Cr(:,:), dimensionless_Fo(:,:), dimensionless_Fi(:,:)
    double precision, allocatable :: dimensionless_Di(:,:), dimensionless_Fc(:,:), dimensionless_D(:,:)

    double precision, allocatable :: ini_y(:), ini_q(:)
    double precision, allocatable :: ini_q_repeat(:,:), ini_E(:,:), ini_F(:,:)

    integer, allocatable :: Q_sk_tableEntry(:,:), noLatFlow(:), noQSKtable(:)
    double precision, allocatable :: eachQSKtableNodeRange(:,:,:), Q_sk_Table(:,:,:,:)

    double precision, allocatable :: lowerLimitCount(:), higherLimitCount(:), volRemain(:,:)

    character(len=128), allocatable :: downstream_path(:), xSection_path(:), manning_strickler_path(:), upstream_path(:),dx_path(:)
    character(len=128), allocatable :: QSKtablePath(:), lateralFlow_path(:), lateralFlow_path2(:)
    character(len=128), allocatable :: bankLocation_path(:)
    double precision, allocatable :: leftBank(:,:), rightBank(:,:), skLeft(:,:), skMain(:,:), skRight(:,:)

    integer, allocatable :: currentROutingDiffusive(:), notSwitchRouting(:)
    double precision :: minDx, maxCelerity,maxCelDx

    integer, allocatable :: currentRoutingNormal(:,:), routingNotChanged(:,:)


end module arrays_module
