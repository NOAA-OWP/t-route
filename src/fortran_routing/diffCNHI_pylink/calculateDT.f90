!Written  by Md Nazmul Azim Beg*, Ehab A Meselhe**
!*:     Postdoctoral Fellow, Tulane River and Coastal Center,
!       Department of River-Coastal Science and Engineering,
!       Tulane University
!**:    Professor, Department of River-Coastal Science and Engineering,
!       Tulane University

!Modified by Dong Ha Kim, NOAA's Office of Water Prediction, National Water Center
module simtime
    use var
    use nrtype
    implicit none

contains
    subroutine cal_dt

        implicit none
        integer(KIND=i4b) :: a, b
        !! initialTime is in hours
        !! tfin is in hours
        !! time is in minutes
        !! dtini is in seconds
        !! saveInterval is in seconds
        dtini = cfl*dx_mn/cel_mx

        a = floor( (tc - t0*60.) /( saveInterval/60. ))
        b = floor( (tc-t0*60.) +dtini/60.)/( saveInterval/60. )

        if (b .gt. a) then
            dtini = (a+1) * ( saveInterval ) - (tc-t0*60.)*60.
        end if

        if ( tc+dtini/60. .gt. tfin*60. ) dtini =  (tfin*60.-tc)*60.

    end subroutine cal_dt
end module simtime
