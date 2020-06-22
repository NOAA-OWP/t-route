module module_reservoir
  
    use precis
    implicit none

contains

subroutine run_reservoir_interface(previous_timestep_inflow, &
    inflow, &
    lateral_inflow, water_elevation, outflow, routing_period, &
    dynamic_reservoir_type)

    implicit none

    !import reservoir
    !class(reservoir), intent(inout) :: this
    real(prec), intent(in) :: dt 
    real(prec), intent(in)    :: previous_timestep_inflow ! cubic meters per second (cms)
    real(prec), intent(in)    :: inflow                   ! cubic meters per second (cms)
    real(prec), intent(in)    :: lateral_inflow           ! cubic meters per second (cms)
    real(prec), intent(inout) :: water_elevation          ! meters AMSL
    real(prec), intent(out)   :: outflow                  ! cubic meters per second (cms)
    real(prec), intent(in)    :: routing_period           ! seconds
    integer, intent(out):: dynamic_reservoir_type   ! dynamic reservoir type sent to lake out files

end subroutine run_reservoir_interface

end module module_reservoir
