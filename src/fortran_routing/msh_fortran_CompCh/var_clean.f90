module var

    implicit none
    save

    integer :: ncomp, ntim, mxncomp
    integer :: ots, option_dsbc
    !* variables for branching channels
    integer :: nlinks, NAnum, mxnbrch, nends
    real :: cfl, f, thes, phi, theta, thetas, thesinv, alfa2, alfa4
    real :: yy, skk, qq, yds
    real :: dtini, dxini, tfin

end module var
