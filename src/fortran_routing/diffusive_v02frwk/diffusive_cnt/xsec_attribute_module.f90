module xsec_attribute_module

    implicit none
    save
    double precision, allocatable :: xsec_tab(:, :, :, :)

    contains
        ! Allocate storage for all of the arrays in this module based on the number
        ! of time steps and spatial points
        subroutine setup_xsec_attribute_module(elements, num_points, num_reaches)
            implicit none
            integer, intent(in) :: elements, num_points, num_reaches
            allocate(xsec_tab(11, elements, num_points, num_reaches))
        end subroutine
end module xsec_attribute_module
