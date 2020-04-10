module arrays

    implicit none
    save

    real, allocatable :: ci1(:), ci2(:)
    real, allocatable :: aso(:), f1(:), f2(:), depth(:)
    real, allocatable :: g11inv(:), g12inv(:), g21inv(:), g22inv(:)
    real, allocatable :: b11(:), b12(:), b21(:), b22(:)
    real, allocatable :: eps2(:), eps4(:), d1(:), d2(:), u(:), c(:)
    real, allocatable :: co(:), gso(:), dbdx(:)
    real, allocatable :: dt(:)
    integer, allocatable :: ityp(:)
    !**arrays for branching channel application
    integer, allocatable :: ndep(:), uslinks(:,:), dslink(:), instrdflag(:,:), nx1(:)
    real, allocatable :: area(:,:), areafnal(:,:,:), y(:, :, :), q(:, :, :), qlat(:,:,:), bo(:, :), traps(:,:)
    real, allocatable :: areap(:, :), qp(:, :), z(:, :), sk(:, :), dx(:, :)
    real, allocatable :: dqp(:, :), dap(:, :), dqc(:, :), dac(:, :)
    real, allocatable :: Tw(:,:), TwCC(:,:), dbdxCCm(:), dbdxCCf(:), skCC1(:,:), skCC2(:,:)
    !integer, allocatable :: chshp(:)

contains

    ! Allocate storage for all of the arrays in this module based on the number
    ! of time steps and spatial points
    subroutine setup_arrays(num_time, mxncomp, nlinks, mxnbrch, nends)

        implicit none
        ! Input
        integer, intent(in) :: num_time, mxncomp, nlinks, mxnbrch, nends

        allocate(area(mxncomp, nlinks), areafnal(num_time,mxncomp,nlinks))
        allocate(y(num_time, mxncomp, nlinks))
        allocate(q(num_time, mxncomp, nlinks))
        allocate(qlat(num_time, mxncomp, nlinks))
        allocate(bo(mxncomp, nlinks), traps(mxncomp,nlinks))

        allocate(Tw(mxncomp, nlinks), TwCC(mxncomp, nlinks), dbdxCCm(mxncomp), dbdxCCf(mxncomp))
        allocate(skCC1(mxncomp, nlinks), skCC2(mxncomp, nlinks))
        !allocate(chshp(mxncomp))

        allocate(areap(mxncomp, nlinks))
        allocate(qp(mxncomp, nlinks))
        allocate(z(mxncomp, nlinks))
        allocate(dqp(mxncomp, nlinks))
        allocate(dqc(mxncomp, nlinks))
        allocate(dap(mxncomp, nlinks))
        allocate(dac(mxncomp, nlinks))
        allocate(ci1(mxncomp))
        allocate(ci2(mxncomp))
        allocate(aso(mxncomp))
        allocate(depth(mxncomp))
        allocate(f1(mxncomp))
        allocate(f2(mxncomp))
        allocate(g11inv(mxncomp))
        allocate(g12inv(mxncomp))
        allocate(g21inv(mxncomp))
        allocate(g22inv(mxncomp))
        allocate(b11(mxncomp))
        allocate(b12(mxncomp))
        allocate(b21(mxncomp))
        allocate(b22(mxncomp))
        allocate(eps2(mxncomp))
        allocate(eps4(mxncomp))
        allocate(d1(mxncomp))
        allocate(d2(mxncomp))
        allocate(u(mxncomp))
        allocate(c(mxncomp))
        allocate(sk(mxncomp, nlinks))
        allocate(co(mxncomp))
        allocate(gso(mxncomp))
        allocate(dbdx(mxncomp))
        allocate(dt(mxncomp))
        allocate(ityp(mxncomp))
        allocate(dx(mxncomp, nlinks))
        allocate(ndep(nlinks), uslinks(mxnbrch,nlinks), dslink(nlinks), nx1(nlinks))
        allocate(instrdflag(nlinks,nends))

    end subroutine setup_arrays

end module arrays
