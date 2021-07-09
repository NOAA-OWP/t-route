module diffusive_interface

use, intrinsic :: iso_c_binding
use diffusive, only: diffnw

implicit none
contains
subroutine c_diffnw(dtini_g, t0_g, tfin_g, saveinterval_g, saveinterval_ev_g, dt_ql_g, dt_ub_g,&
        dt_db_g, nts_ql_g, nts_ub_g, nts_db_g, mxncomp_g, nrch_g, z_ar_g, bo_ar_g, traps_ar_g,&
        tw_ar_g, twcc_ar_g, mann_ar_g, manncc_ar_g, so_ar_g, dx_ar_g, iniq, nhincr_m_g, nhincr_f_g,&
        ufhlt_m_g,  ufqlt_m_g, ufhlt_f_g, ufqlt_f_g, frnw_col, frnw_g, qlat_g, ubcd_g, dbcd_g,&
        cfl_g, theta_g, tzeq_flag_g, y_opt_g, so_llm_g, ntss_ev_g, q_ev_g, elv_ev_g) bind(c)

    real(c_double), intent(in) :: dtini_g, t0_g, tfin_g, saveinterval_g, saveinterval_ev_g
    real(c_double), intent(in) :: dt_ql_g, dt_ub_g, dt_db_g
    integer(c_int), intent(in) :: nts_ql_g, nts_ub_g, nts_db_g
    integer(c_int), intent(in) :: mxncomp_g, nrch_g
    real(c_double), dimension(mxncomp_g, nrch_g), intent(in) :: z_ar_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g
    real(c_double), dimension(mxncomp_g, nrch_g), intent(in) :: mann_ar_g, manncc_ar_g
    real(c_double), dimension(mxncomp_g, nrch_g), intent(inout) :: so_ar_g
    real(c_double), dimension(mxncomp_g, nrch_g), intent(in) :: dx_ar_g, iniq
    integer(c_int), intent(in) :: nhincr_m_g, nhincr_f_g
    real(c_double), dimension(mxncomp_g, nrch_g, nhincr_m_g), intent(in) :: ufhlt_m_g,  ufqlt_m_g
    real(c_double), dimension(mxncomp_g, nrch_g, nhincr_f_g), intent(in) :: ufhlt_f_g, ufqlt_f_g
    integer(c_int), intent(in) :: frnw_col
    real(c_double), dimension(nrch_g, frnw_col), intent(in) :: frnw_g
    real(c_double), dimension(nts_ql_g, mxncomp_g, nrch_g), intent(in) :: qlat_g
    real(c_double), dimension(nts_ub_g, nrch_g), intent(in) :: ubcd_g
    real(c_double), dimension(nts_db_g), intent(in) :: dbcd_g
    real(c_double), intent(in) :: cfl_g, theta_g
    integer(c_int), intent(in) :: tzeq_flag_g
    integer(c_int), intent(in) :: y_opt_g
    real(c_double), intent(in) :: so_llm_g
    integer(c_int), intent(in) :: ntss_ev_g
    real(c_double), dimension(ntss_ev_g, mxncomp_g, nrch_g), intent(out) :: q_ev_g, elv_ev_g
    
    call diffnw(dtini_g, t0_g, tfin_g, saveinterval_g, saveinterval_ev_g, dt_ql_g, dt_ub_g,&
        dt_db_g, nts_ql_g, nts_ub_g, nts_db_g, mxncomp_g, nrch_g, z_ar_g, bo_ar_g, traps_ar_g,&
        tw_ar_g, twcc_ar_g, mann_ar_g, manncc_ar_g, so_ar_g, dx_ar_g, iniq, nhincr_m_g, nhincr_f_g,&
        ufhlt_m_g,  ufqlt_m_g, ufhlt_f_g, ufqlt_f_g, frnw_col, frnw_g, qlat_g, ubcd_g, dbcd_g,&
        cfl_g, theta_g, tzeq_flag_g, y_opt_g, so_llm_g, ntss_ev_g, q_ev_g, elv_ev_g)
    
end subroutine c_diffnw
end module diffusive_interface
