module var
    !* This module is for F2pying MSH model on multiple segments of a given single reach at a given single time.
    !* p.F<->P,RM3

    implicit none
!    save



    !integer :: ncomp, ntim, mxncomp
    !integer :: ots, option_dsbc
    !* variables for branching channels
    !integer :: nlinks, NAnum, mxnbrch, nends
    !real :: cfl, f, thes, phi, theta, thetas, thesinv, alfa2, alfa4
    !real :: yy, skk, qq, yds
    !real :: dtini, dxini, tfin
    real, parameter :: grav = 9.81
    real, parameter :: TOLERANCE = 1e-8
    !*-------------------------------------------------------------
    !* MSH coefficients
    real :: cfl, f, thes, phi, theta, thetas, thesinv, alfa2, alfa4
    integer :: ots
    integer, allocatable :: ityp(:)
    !*-------------------------------------------------------------

    !integer :: j, n
    integer :: ncompy !  ntim, nlinks, mxncomp
    !* variables for branching channels
    !integer :: NAnum, mxnbrch, nends
    !real :: yy,  qq, yds, yini !, skk
    real :: dtini !, dxini, tfin
    integer :: eqn_A

    !integer, allocatable :: ndep(:), uslinks(:,:), dslink(:), instrdflag(:,:), nx1(:)
    !real, allocatable, dimension(:,:,:) :: y, q, qlat

    !* channel data
    real, allocatable, dimension(:) :: dx, z, bo, traps, sk, chbtslp
    real, allocatable, dimension(:) :: Tw, TwCC, skCC1, skCC2
    !* lateral inflow
    real, allocatable, dimension(:) :: qlat
    !* time step
    real, allocatable, dimension(:) :: dt
    !* MSH computation variables
    real, allocatable, dimension(:) :: y, q, area , dap, dqp, dac, dqc, areap, qp
    real, allocatable, dimension(:) :: ynp1, qnp1 ! ynm1 !* y(n+1), q(n+1), y(n-1)
    real :: ynm1 !*for upstream boundary condition
    !* boundary conditions
    integer :: usbdflag, dsbdflag !i_usbd,
    real :: msrq, msrqnp1
    real :: msrh, msrhnp1

    !* for test
!    real, allocatable, dimension(:) :: cop, ci1p, gsop, depthp, b11p, b12p, b21p, b22p
!    real, allocatable, dimension(:) :: g11invp, g12invp, g21invp, g22invp, f1p, f2p, d1p, d2p
!    real, allocatable, dimension(:) :: up, cp, a11p, a12p, a21p, a22p, st11p, st12p ,st21p, st22p
!    real, allocatable, dimension(:,:) :: courp


    !real :: sbtoolrst  !*represent a result of any procedures of sub subtools
                       !*that want to make results accessible globally.
    !real :: estqn !* estimated q at time n (used for initial condition)
    !real :: estdepthn !* estimated depth at time n (used for initial condition)
    !integer :: usbd_flag, dsbd_flag
    !integer :: nseg !* the number of segments of a given reach
    !real :: eqnp1 !* =q(n+1,1,j)
    !real :: eynp1 !* =y(n+1,ncomp,j)
    !real, allocatable, dimension(:) :: ynp1, qnp1 !*=y(n+1,i,j) and q(n+1,i,j), respectively.

end module var




