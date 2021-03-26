module mdiffv2_interface

use, intrinsic :: iso_c_binding
use mdiffv2, only: diffnw

implicit none
contains
subroutine c_diffnw(mxncomp_g, nrch_g, z_ar_g, ntss_ev_g, q_ev_g, elv_ev_g) bind(c)

    integer, intent(in) :: mxncomp_g, nrch_g, ntss_ev_g
    double precision, dimension(mxncomp_g, nrch_g), intent(in) :: z_ar_g
    double precision, dimension(ntss_ev_g, mxncomp_g, nrch_g), intent(out) :: q_ev_g, elv_ev_g

    call diffnw(mxncomp_g, nrch_g, z_ar_g, ntss_ev_g, q_ev_g, elv_ev_g)
    
end subroutine c_diffnw
end module mdiffv2_interface