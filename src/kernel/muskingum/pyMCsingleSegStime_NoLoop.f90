module muskingcunge_interface

use, intrinsic :: iso_c_binding, only: c_float
use muskingcunge_module, only: muskingcungenwm, reachcompute

implicit none
contains
subroutine c_muskingcungenwm(dt, qup, quc, qdp, ql, dx, bw, tw, twcc,&
    n, ncc, cs, s0, velp, depthp, qdc, velc, depthc, ck, cn, X) bind(c)

    real(c_float), intent(in) :: dt
    real(c_float), intent(in) :: qup, quc, qdp, ql
    real(c_float), intent(in) :: dx, bw, tw, twcc, n, ncc, cs, s0
    real(c_float), intent(in) :: velp, depthp
    real(c_float), intent(out) :: qdc, velc, depthc
    real(c_float), intent(out) :: ck, cn, X

    call muskingcungenwm(dt, qup, quc, qdp, ql, dx, bw, tw, twcc,&
    n, ncc, cs, s0, velp, depthp, qdc, velc, depthc, ck, cn, X)
    !print*, "fortran c_bind", depthc
    
end subroutine c_muskingcungenwm

subroutine c_reachcompute(dt, nseg, nts, qup_top, quc_top, qdp_rch, ql_rch, dx_rch,& 
                    bw_rch, tw_rch, twcc_rch,n_rch, ncc_rch, cs_rch, s0_rch,&
                    velp_rch, depthp_rch, qdc_rch, velc_rch, depthc_rch, ck_rch,&
                    cn_rch, X_rch) bind(c)
                    
    integer, intent(in) :: nseg, nts
    real(c_float), intent(in) :: dt
    real(c_float), dimension(nts), intent(in) :: qup_top, quc_top
    real(c_float), dimension(nseg, nts), intent(in) :: ql_rch
    real(c_float), dimension(nseg), intent(in) :: qdp_rch
    real(c_float), dimension(nseg), intent(in) :: dx_rch, bw_rch, tw_rch, twcc_rch
    real(c_float), dimension(nseg), intent(in) :: n_rch, ncc_rch, cs_rch, s0_rch
    real(c_float), dimension(nseg), intent(in) :: velp_rch
    real(c_float), dimension(nseg), intent(in) :: depthp_rch
    real(c_float), dimension(nseg, nts), intent(out) :: qdc_rch, velc_rch, depthc_rch
    real(c_float), dimension(nseg, nts), intent(out) :: ck_rch, cn_rch, X_rch
                        
    call reachcompute(dt, nseg, nts, qup_top, quc_top, qdp_rch, ql_rch, dx_rch,& 
                      bw_rch, tw_rch, twcc_rch,n_rch, ncc_rch, cs_rch, s0_rch,&
                      velp_rch, depthp_rch, qdc_rch, velc_rch, depthc_rch, ck_rch,&
                      cn_rch, X_rch)
                        
end subroutine c_reachcompute

end module muskingcunge_interface
