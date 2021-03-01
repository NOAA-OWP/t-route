module gfunc_module
implicit none
    double precision :: mt
contains
subroutine gfunc(x, n, m, a, b, c)
    double precision, intent(in) :: x
    integer, intent(in) :: n, m
    double precision, dimension(n), intent(in) :: a
    double precision, dimension(m), intent(in) :: b
    double precision, dimension(n, m), intent(out) :: c
    integer :: i, j

    do j=1,m
        do i=1,n
             c(i,j) = exp(-x * (a(i)**2 + b(j)**2))
        end do
    end do
    mt=2.0
    call gfunc2(n,m,c)

end subroutine

subroutine gfunc2(n, m, c)
    integer, intent(in) :: n, m
    double precision, dimension(n,m), intent(inout) :: c
    integer :: i,j

    do j=1,m
        do i=1,n
            c(i,j)= mt*c(i,j)
        enddo
    enddo  
end subroutine gfunc2
end module
