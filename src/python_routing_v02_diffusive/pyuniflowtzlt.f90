module uniflowtzlt_interface
use iso_c_binding, only: c_double, c_int
use flowlt, only: uniflowtzlt
implicit none
contains
subroutine c_uniflowtzlt(mxncomp_g, nrch_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g, &
                        mann_ar_g, manncc_ar_g, so_ar_g, &
                        nhincr_m_g, nhincr_f_g, &
                        ufhlt_m_g, ufqlt_m_g, ufhlt_f_g, ufqlt_f_g, &
                        frnw_col, dfrnw_g, timesdepth_g) bind(c)
    integer(c_int), intent(in) :: mxncomp_g, nrch_g
    integer(c_int), intent(in) :: nhincr_m_g, nhincr_f_g, frnw_col
    real(c_double), dimension(mxncomp_g, nrch_g), intent(in) :: bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g
    real(c_double), dimension(mxncomp_g, nrch_g), intent(in) :: mann_ar_g, manncc_ar_g, so_ar_g
    real(c_double), dimension(nrch_g, frnw_col), intent(in) :: dfrnw_g 
    real(c_double), intent(in) :: timesdepth_g
    real(c_double), dimension(mxncomp_g, nrch_g, nhincr_m_g), intent(out) :: ufhlt_m_g,  ufqlt_m_g
    real(c_double), dimension(mxncomp_g, nrch_g, nhincr_f_g), intent(out) :: ufhlt_f_g, ufqlt_f_g
    call uniflowtzlt(mxncomp_g, nrch_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g, &
                     mann_ar_g, manncc_ar_g, so_ar_g, &
                     nhincr_m_g, nhincr_f_g, &
	     	     ufhlt_m_g, ufqlt_m_g, ufhlt_f_g, ufqlt_f_g, &
                     frnw_col, dfrnw_g, timesdepth_g)
end subroutine
end module
