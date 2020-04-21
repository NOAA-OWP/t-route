module var
    implicit none
    save

    real :: dt, qup, quc, qdp, qdc, ql, z, vel, depth
    real :: bfd, WPC, AREAC, C1, C2, C3, C4

    integer :: ntim
    !integer :: ncomp0, ncomp, uslinkID, linkID !,nlinks
    integer :: ncomp, linkID !,nlinks

    real,allocatable,dimension(:) :: dx, bw, tw, twcc, ncc, cs, so, n
    real,allocatable,dimension(:) :: qlat !real,allocatable,dimension(:,:) :: qlat
    real,allocatable,dimension(:,:) :: vela, deptha
    real,allocatable,dimension(:,:,:) :: qd

end module var



