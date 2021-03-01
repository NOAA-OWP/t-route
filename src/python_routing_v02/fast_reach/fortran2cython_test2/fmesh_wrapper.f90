module fmesh_wrapper

  use iso_c_binding, only: c_double, c_int
  use fmesh, only: meshexp

  implicit none

contains

  subroutine c_meshexp(rmin, rmax, a, N, mesh) bind(c)
    real(c_double), intent(in) :: rmin
    real(c_double), intent(in) :: rmax
    real(c_double), intent(in) :: a
    integer(c_int), intent(in) :: N
    real(c_double), intent(out) :: mesh(N+1)
    mesh = meshexp(rmin, rmax, a, N)
  end subroutine c_meshexp

  ! wrap more functions here
  ! ...

end module fmesh_wrapper
