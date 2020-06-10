module var
    implicit none
    save

    integer, parameter :: dp = kind(1.d0)
    real(dp) :: dt, dx, qup, quc, qdp, qdc, ql
    real(dp) :: Bw, Tw, TwCC, nCC, Cs, So, n, z, vel, depth
    real(dp) :: velp_chk, depthp_chk
    real(dp) :: bfd, WPC, AREAC, C1, C2, C3, C4
    integer :: ntim
    integer :: ncomp0, ncomp, iseg, uslinkflag
    real(dp) :: Qus_prev
    real(dp),allocatable,dimension(:,:) :: vela, deptha
    real(dp),allocatable,dimension(:,:) :: Qd


end module var



