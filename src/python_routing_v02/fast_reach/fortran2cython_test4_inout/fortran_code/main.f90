module mdiffv2
contains
    
    subroutine diffnw(mxncomp_g, nrch_g, z_ar_g, ntss_ev_g, q_ev_g, elv_ev_g)
    
        implicit none
        integer, intent(in) :: mxncomp_g, nrch_g, ntss_ev_g
        double precision, dimension(mxncomp_g, nrch_g), intent(in) :: z_ar_g
        double precision, dimension(ntss_ev_g, mxncomp_g, nrch_g), intent(out) :: q_ev_g, elv_ev_g
        
        integer :: plane
        integer :: row
        integer :: col
        
        do plane = 1, ntss_ev_g
            do row = 1, mxncomp_g
                do col = 1, nrch_g
                    q_ev_g(plane,row,col) = z_ar_g(row,col)
                    elv_ev_g(plane,row,col) = z_ar_g(row,col) * 2.0
            end do
          end do
        end do

    end subroutine diffnw
end module mdiffv2