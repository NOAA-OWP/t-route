module muskingcunge_interface

use iso_c_binding, only: c_float
use muskingcunge_module, only: muskingcungenwm

implicit none
contains

subroutine c_muskingcungenwm(dt, qup, quc, qdp, ql, dx, bw, tw, twcc,&
    n, ncc, cs, s0, velp, depthp, qdc, velc, depthc) bind(c)

    real(c_float), intent(in) :: dt
    real(c_float), intent(in) :: qup, quc, qdp, ql
    real(c_float), intent(in) :: dx, bw, tw, twcc, n, ncc, cs, s0
    real(c_float), intent(in) :: velp, depthp
    real(c_float), intent(out) :: qdc, velc, depthc

    call muskingcungenwm(dt, qup, quc, qdp, ql, dx, bw, tw, twcc,&
    n, ncc, cs, s0, velp, depthp, qdc, velc, depthc)
end subroutine

end module