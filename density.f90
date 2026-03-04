module get_density
    implicit none
CONTAINS
    subroutine get_den(n, nghost, x, m, h, rho)
        real, dimension(:), intent(in) :: x, m, h
        integer, intent(in) :: n, nghost
        real, dimension(:), intent(inout) :: rho
        integer :: i, j
        real :: dist, q
        real, PARAMETER :: sigma = 2.0/3.0

        do i=1, n
            rho(i) = 0.0
            do j=1,nghost+n
                dist = abs(x(i)-x(j))
                q = dist/h(i)
                rho(i) = rho(i) + m(j) * sigma/h(i) * W_calc(q)
            end do
        end do

    end subroutine get_den

    pure function W_calc(q) result(w)
        real, intent(IN) :: q
        real :: w

        if (q >= 0.0 .AND. q <1.0) THEN
            w = 0.25 * (2.0-q)**3 - (1.0-q)**3
        elseif (q >=1.0 .AND. q<2.0) THEN
            w = 0.25 * (2.0-q)**3
        else
            w = 0.0
        end if
    end function W_calc

end module get_density

! Check which dimension sigma wants us to use