module var
    use precis
    implicit none
    save

    real(prec) :: dt, dx, qup, quc, qdp, qdc, ql
    real(prec) :: Bw, Tw, TwCC, nCC, Cs, So, n, z, vel, depth
    real(prec) :: velp_chk, depthp_chk
    real(prec) :: bfd, WPC, AREAC, C1, C2, C3, C4
    integer :: ntim
    integer :: ncomp0, ncomp, iseg, uslinkflag
    real(prec) :: Qus_prev
    real(prec),allocatable,dimension(:,:) :: vela, deptha
    real(prec),allocatable,dimension(:,:) :: Qd


end module var



