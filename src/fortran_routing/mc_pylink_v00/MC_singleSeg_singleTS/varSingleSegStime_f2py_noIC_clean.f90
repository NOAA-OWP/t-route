module var
    implicit none
    save

    integer, parameter :: dp = kind(1.d0)
    real(dp) :: dt, dx, qup, quc, qdp, qdc, ql, Bw, Tw, TwCC, nCC, Cs, So, n, z, vel, depth
    real(dp) :: bfd, WPC, AREAC, C1, C2, C3, C4

end module var



