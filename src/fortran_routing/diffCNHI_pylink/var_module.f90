module var

    use nrtype

    implicit none

    !save

    real(KIND=dp), parameter :: grav = 9.81
    real(KIND=dp), parameter :: TOLERANCE = 1e-8
    !character*4 :: file_num
    !character(len=128) :: xSection_path

    !integer(KIND=i4b) :: ncomp_unif_g
    integer(KIND=i4b) :: nel_g      !* # of potential elevation lines for x-sec attr.table
    integer(KIND=i4b) :: nxsecpt_g  !* # of x-sec x-y data points
    !real(KIND=dp) :: f_g
    real(KIND=dp) :: dtini_g, t0_g, tfin_g,  cfl_g
    real(KIND=dp) :: saveInterval_g, saveInterval_ev_g

    real(KIND=dp), dimension(:), allocatable :: xcs_g, ycs_g
    !real(KIND=dp), allocatable :: xsec_tab(:, :, :)
    real(KIND=dp), allocatable :: xsec_attr_seg_g(:,:) !* first col: type of attribute(dim:3)
                                            !* second col: elev increments at each segment (dim:nel)
    real(KIND=dp), allocatable :: xsec_attr_rch_g(:,:,:,:)   !* first col: type of attribute (dim:3)
                                                !* second col: elev increments (dim:nel)
                                                !* third col: segID of a given reach (dim:ncomp)
    !real(KIND=dp), dimension(:), allocatable :: dx_ar_g !,  sk ! z,
    !real(KIND=dp), dimension(:), allocatable :: celty_g, diffty_g
    !real(KIND=dp), dimension(:), allocatable :: q, qpx_g !* when passed by python, the values are of current time tc,
                                              !* but when being passed to python from fortran, the values
                                              !* are of the next time (= tc + dtini/60 [min])
    !real(KIND=dp), dimension(:), allocatable :: elv    !* Python only received computed values of elv by Fortran, which
                                              !* are of tc + dtini/60 [min]

    real(KIND=dp) :: timesDepth_g  !* multiplier to potentially largest depth

    !* for analytically computing channel properties
    integer(KIND=i4b) :: tzeq_flag_g

    real(KIND=dp):: z_g, bo_g, traps_g, tw_g, twcc_g, So_g, mann_g, manncc_g
    !real(KIND=dp), dimension(:), allocatable :: z_r, bo_r, traps_r, tw_r, twcc_r, So_r
    !real(KIND=dp), dimension(:), allocatable :: mann_r, manncc_r
    !*for subroutine calculateDT
    real(KIND=dp) :: dx_mn_g !, cel_mx
    !real(KIND=dp), allocatable :: qlatj_g(:)
    !real(KIND=dp) :: qlatf
    !real(KIND=dp), allocatable :: cel_av_g(:)
    real(KIND=dp) :: theta_g
    integer(KIND=i4b) :: y_opt_g      !* 1:normal depth elevation using either bisection or lookup table
                                    !* 2: diffusive wave elevation

    !* added variables for Fortran for all project ^^!!
    integer(KIND=i4b) :: mxncomp_g, nrch_g
    integer(KIND=i4b) :: ntss_g, nts_ql_g, nts_ub_g, nts_db_g, ntss_ev_g


    integer(KIND=i4b) :: nhincr_m_g, nhincr_f_g   !*the number of depth increment of main ch. and floodplain.
    real(KIND=dp) :: dt_ql_g, dt_ub_g, dt_db_g
    !real(KIND=dp), allocatable :: xsec_attr_rch(:,:,:,:)
    !real(KIND=dp), dimension(:,:), allocatable :: ufhlt_mr,  ufqlt_mr, ufhlt_fr, ufqlt_fr
    real(KIND=dp), dimension(:,:,:), allocatable :: ufhlt_m_g,  ufqlt_m_g, ufhlt_f_g, ufqlt_f_g
    integer(KIND=i4b), dimension(:,:), allocatable :: frnw_g
    real(KIND=dp), dimension(:,:,:), allocatable :: qlat_g
    real(KIND=dp), dimension(:,:), allocatable :: ubcd_g
    real(KIND=dp), dimension(:), allocatable :: dbcd_g

    real(KIND=dp), dimension(:,:), allocatable :: z_ar_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g, So_ar_g
    real(KIND=dp), dimension(:,:), allocatable :: mann_ar_g, manncc_ar_g, dx_ar_g

    real(KIND=dp), dimension(:,:,:), allocatable :: q_ev_g, elv_ev_g !,q, elv_g
    real(KIND=dp), dimension(:,:), allocatable :: celty_g, diffty_g, q_g, elv_g, qpx_g
    real(KIND=dp), dimension(:), allocatable :: cel_av_g
    real(KIND=dp), dimension(:), allocatable :: qlatj_g
    real(KIND=dp), dimension(:), allocatable :: eei_g, ffi_g, exi_g, fxi_g
    !* variable for the algorithm with modified segment lengths
    real(KIND=dp), dimension(:), allocatable :: celty_m_g, diffty_m_g, q_m_g, qpx_m_g
    real(KIND=dp), dimension(:), allocatable :: dx_ar_m_g
    real(KIND=dp), dimension(:), allocatable :: qlatj_m_g
    real(KIND=dp), dimension(:), allocatable :: eei_m_g, ffi_m_g, exi_m_g, fxi_m_g

    real(KIND=dp) :: so_llm_g
    integer(KIND=i4b) :: nl_ubcd_ar_g


end module var
