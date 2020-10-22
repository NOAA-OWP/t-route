module muskingcunge_interface

use, intrinsic :: iso_c_binding, only: c_float
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
<<<<<<< HEAD

    call muskingcungenwm(dt, qup, quc, qdp, ql, dx, bw, tw, twcc,&
    n, ncc, cs, s0, velp, depthp, qdc, velc, depthc)
end subroutine c_muskingcungenwm
end module muskingcunge_interface
=======
    real(c_float) :: ck, cn, X
    !TODO: Incorporate ck, cn, X into v02 output;
    ! these are currently dropped silently from
    ! the interface output.

    call muskingcungenwm(dt, qup, quc, qdp, ql, dx, bw, tw, twcc,&
    n, ncc, cs, s0, velp, depthp, qdc, velc, depthc, ck, cn, X)
    
end subroutine c_muskingcungenwm
end module muskingcunge_interface
>>>>>>> upstream/master
