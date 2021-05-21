module dimlessN
    use var
    implicit none

contains
    subroutine calc_dimensionless_numbers(j)
        implicit none

        integer, intent(in) :: j
        integer(kind=4) :: i
        real(kind=4) :: wl_us, depth_us, q_us, v_us, pere_us, r_us, sk_us, ch_us
        real(kind=4) :: wl_ds, depth_ds, q_ds, v_ds, pere_ds, r_ds, sk_ds, ch_ds
        real(kind=4) :: ch_star_avg, channel_length, avg_celerity, avg_velocity, avg_depth
        real(kind=4) :: maxValue, dimlessWaveLength
        maxValue = 1e7

        do i=1, ncomp-1

            avg_celerity = (celerity(i,j) + celerity(i+1,j)) / 2.0    ! 'celerity2' is calculated celerity. 'celerity' is spatially averaged

            wl_us = newY(i,j)
            ! bo(i) and pere(i) has data for the latest river reach only
            depth_us = newArea(i,j) / bo(i,j)
            q_us = newQ(i,j)
            v_us = abs( newQ(i,j) / newArea(i,j) )
            pere_us = pere(i,j)
            r_us = newArea(i,j) / pere(i,j)
            sk_us = sk(i,j)
            ch_us = sk(i,j) * r_us ** (1./6.)


            wl_ds = newY(i+1,j)
            depth_ds = newArea(i+1,j) / bo(i+1,j)
            q_ds = newQ(i+1,j)
            v_ds = abs( newQ(i+1,j) / newArea(i+1,j) )
            pere_ds = pere(i+1,j)
            r_ds = newArea(i+1,j) / pere(i+1,j)
            sk_ds = sk(i+1,j)
            ch_ds = sk(i+1,j) * r_ds ** (1./6.)


            ch_star_avg = ((ch_us + ch_ds) / 2.)  / sqrt( grav ) !! CORRECTED
            channel_length = dx(i,j)

            dimlessWaveLength= 4000. !! new

            avg_velocity = (v_us + v_ds) / 2.
            avg_depth = (depth_us + depth_ds) / 2.

            dimensionless_Cr(i,j) = abs(avg_velocity / avg_celerity)
            if (dimensionless_Cr(i,j) .gt. maxValue) dimensionless_Cr(i,j) = maxValue

            dimensionless_Fo(i,j) = avg_velocity / sqrt(grav * avg_depth)
            if (dimensionless_Fo(i,j) .gt. maxValue) dimensionless_Fo(i,j) = maxValue

            !dimensionless_Fi(i) = 2*dimensionless_Cr(i) / (ch_star_avg ** 2.) * (channel_length / avg_depth) !! CORRECTED
            dimensionless_Fi(i,j) = 2*dimensionless_Cr(i,j) / (ch_star_avg ** 2.) * dimlessWaveLength !(channel_length / avg_depth) !! CORRECTED
            if (dimensionless_Fi(i,j) .gt. maxValue) dimensionless_Fi(i,j) = maxValue

            dimensionless_Fc(i,j) = dimensionless_Cr(i,j) * dimensionless_Fi(i,j)
            if (dimensionless_Fc(i,j) .gt. maxValue) dimensionless_Fc(i,j) = maxValue

            dimensionless_Di(i,j) = (dimensionless_Cr(i,j) / dimensionless_Fo(i,j)) ** 2. !! CORRECTED
            if (dimensionless_Di(i,j) .gt. maxValue) dimensionless_Di(i,j) = maxValue

            dimensionless_D(i,j)  = dimensionless_Di(i,j) / dimensionless_Fc(i,j)
            if (dimensionless_D(i,j) .gt. maxValue) dimensionless_D(i,j) = maxValue
        end do
    end subroutine calc_dimensionless_numbers
endmodule dimlessN
