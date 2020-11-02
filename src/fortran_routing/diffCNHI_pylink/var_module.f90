!Written  by Md Nazmul Azim Beg*, Ehab A Meselhe**
!*:     Postdoctoral Fellow, Tulane River and Coastal Center,
!       Department of River-Coastal Science and Engineering,
!       Tulane University
!**:    Professor, Department of River-Coastal Science and Engineering,
!       Tulane University

!Modified by Dong Ha Kim, NOAA's Office of Water Prediction, National Water Center
module var
    use nrtype
    implicit none
    save

    real(KIND=dp), parameter :: grav = 9.81
    real(KIND=dp), parameter :: TOLERANCE = 1e-8

    integer(KIND=i4b) :: ncomp
    integer(KIND=i4b) :: nel      !* # of potential elevation lines for x-sec attr.table
    integer(KIND=i4b) :: nxsecpt  !* # of x-sec x-y data points
    real(KIND=dp) :: cfl, f
    real(KIND=dp) :: dtini, t0, tc, saveInterval, tfin

    real(KIND=dp), dimension(:), allocatable :: xcs, ycs
    real(KIND=dp), allocatable :: xsec_attr_seg(:,:) !* first col: type of attribute(dim:3)
                                            !* second col: elev increments at each segment (dim:nel)
    real(KIND=dp), allocatable :: xsec_attr_rch(:,:,:)   !* first col: type of attribute (dim:3)
                                                !* second col: elev increments (dim:nel)
                                                !* third col: segID of a given reach (dim:ncomp)
    real(KIND=dp), dimension(:), allocatable :: dx !,  sk ! z,
    real(KIND=dp), dimension(:), allocatable :: celty, diffty
    real(KIND=dp), dimension(:), allocatable :: q, qpx !* when passed by python, the values are of current time tc,
                                              !* but when being passed to python from fortran, the values
                                              !* are of the next time (= tc + dtini/60 [min])
    real(KIND=dp), dimension(:), allocatable :: elv    !* Python only received computed values of elv by Fortran, which
                                              !* are of tc + dtini/60 [min]

    real(KIND=dp) :: timesDepth  !* multiplier to potentially largest depth
    real(KIND=dp) :: eqnp1, eynp1 !* estimated q and y(=elevation) at n+1(= t + dtini/60), where t=tc

    integer(KIND=i4b) :: tzeq_flag
    real(KIND=dp):: z0, bo0, traps0, tw0, twcc0, So0, mann, manncc
    real(KIND=dp), dimension(:), allocatable :: z_ar, bo_ar, traps_ar, tw_ar, twcc_ar, So_ar
    real(KIND=dp), dimension(:), allocatable :: mann_ar, manncc_ar
    !*for subroutine calculateDT
    real(KIND=dp) :: dx_mn, cel_mx
    real(KIND=dp), allocatable :: qlatj(:)
    real(KIND=dp) :: qlatf
    real(KIND=dp) :: theta
    !* normal depth-discharge lookup table parameters
    integer(KIND=i4b) :: nhincr_m, nhincr_f !*the number of depth increment of main ch. and floodplain.
    real(KIND=dp) :: mtp_hmx        !*multiplier for estimating potential max. depth
    integer(KIND=i4b) :: y_opt      !* 1:normal depth elevation using either bisection or lookup table
                                    !* 2: diffusive wave elevation
    real(KIND=dp), dimension(:,:), allocatable :: ufqlt_m, ufhlt_m, ufqlt_f, ufhlt_f



end module var
