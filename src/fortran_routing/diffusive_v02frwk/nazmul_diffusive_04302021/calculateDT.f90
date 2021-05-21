module tstep
    use var
    implicit none

contains
    subroutine calculateDT(initialTime, time, saveInterval, maxAllowCourantNo, tfin, max_C_dx, given_dt)
        implicit none

        real, intent(in) :: initialTime, time, saveInterval, tfin, given_dt
        real, intent(in) :: maxAllowCourantNo, max_C_dx
        integer          :: a, b

        !! initialTime is in hours
        !! tfin is in hours
        !! time is in minutes
        !! dtini is in seconds
        !! saveInterval is in seconds
        dtini = maxAllowCourantNo/max_C_dx

        a = floor( (time-initialTime*60.) /( saveInterval/60. ))            ! units:: time : min;  ! initialTime : hour ! saveInterval : sec
        b = floor( (time-initialTime*60.) +dtini/60.)/( saveInterval/60. )

        if (b .gt. a) then
            dtini = (a+1) * ( saveInterval ) - (time-initialTime*60.)*60.
        end if

        if ( time+dtini/60. .gt. tfin*60. ) dtini =  (tfin*60.-time)*60.
    end subroutine

    subroutine correctDT(initialTime, time, saveInterval, tfin)
        implicit none

        real(kind=4), intent(in) :: initialTime, time, saveInterval, tfin
        integer(kind=4)  :: a, b

        !! initialTime is in hours
        !! tfin is in hours
        !! time is in minutes
        !! dtini is in seconds
        !! saveInterval is in seconds
        a = floor( (time-initialTime*60.) /( saveInterval/60. ))
        b = floor(( (time-initialTime*60.) +dtini/60.)/( saveInterval/60. ))

        if( b .gt. a ) then
            dtini = ( initialTime*60. + (a+1)*saveInterval/60. -  time ) * 60.
        end if
        if ( a*saveInterval / 60. +initialTime*60.+dtini / 60. .gt. tfin*60. ) then
            dtini = ( tfin*60. - time )*60.
        end if

    end subroutine
endmodule tstep
