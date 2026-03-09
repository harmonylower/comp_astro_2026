module calcs
    use kernal, only:cubic_spline, grad_cubic_spline
    implicit none
contains
    subroutine get_den(n, n_ghost, x, m, h, rho)
        real, dimension(:), intent(in) :: x, m, h
        integer, intent(in) :: n, n_ghost
        real, dimension(:), intent(inout) :: rho
        integer :: i, j
        real :: dist, q, w
        real, parameter :: sigma = 2.0/3.0

        do i=1, n
            rho(i) = 0.0
            do j=1,n_ghost+n
                dist = abs(x(i)-x(j))
                q = dist/h(i)
                call cubic_spline(q, w)
                rho(i) = rho(i) + m(j) * sigma/h(i) * w
            end do
        end do

    end subroutine get_den

    subroutine get_accel(n, n_ghost, x, m, h, rho, pre, a)
        integer, intent(in) :: n, n_ghost
        real, dimension(:), intent(in) :: x, m, h, rho, pre
        real, dimension(:), intent(inout) :: a

        real :: dist, q_i, q_j, w_i, w_j, lft, rgt
        integer :: i, j
        real, dimension(n+n_ghost*3) :: q_ab
        real, parameter :: sigma = 2.0/3.0
        q_ab(1:n+n_ghost) = 0.0

        a(1:n) = 0.0 
        do i=1, n
            do j=1, n+n_ghost
                if (i==j) cycle !skip itself
                dist = abs(x(i)-x(j))
                q_i = dist/h(i)
                q_j = dist/h(j)
                call grad_cubic_spline(q_i, w_i)
                call grad_cubic_spline(q_j, w_j)
                
                lft = (pre(i) + q_ab(i)) / rho(i)**2  * sigma/h(i)**2 * w_i * sign(1.0,x(i)-x(j))
                rgt = (pre(j) + q_ab(j)) / rho(j)**2  * sigma/h(j)**2 * w_j * sign(1.0,x(i)-x(j))
                a(i) = a(i) - m(j) * (lft + rgt)
                
            end do
        end do
    end subroutine get_accel

    subroutine get_h(m, rho, n, h)
        use setup, only:hfac
        real, dimension(:), intent(in) :: m, rho
        integer, intent(in) :: n
        real, dimension(:), intent(inout) :: h
        integer :: i
    
        do i=1,n
            h(i) = hfac * (m(i)/rho(i)) ! to the power of 1/dim where dim=1
        end do
    end subroutine get_h

end module calcs

!calculating aeleration before ghosts?