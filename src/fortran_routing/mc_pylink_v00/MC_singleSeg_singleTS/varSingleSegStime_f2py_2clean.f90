module var
    implicit none
    save

    real*8 :: dt, dx, qup, quc, qdp, qdc, ql
    real*8 :: Bw, Tw, TwCC, nCC, Cs, So, n, z, vel, depth
    real*8 :: velp_chk, depthp_chk
    real*8 :: bfd, WPC, AREAC, C1, C2, C3, C4
    integer :: ntim
    integer :: ncomp0, ncomp, iseg, uslinkflag
    real*8 :: Qus_prev
    real*8,allocatable,dimension(:,:) :: vela, deptha
    real*8,allocatable,dimension(:,:) :: Qd


end module var



