module calcs
    use kernal, only:cubic_spline, grad_cubic_spline
    use setup, only: fix_choice
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

    subroutine get_accel(n, n_ghost, x, v, m, h, rho, pre, a, cs)
        integer, intent(in) :: n, n_ghost
        real, dimension(:), intent(in) :: x, m, h, rho, pre, v, cs
        real, dimension(:), intent(inout) :: a

        real :: dist, q_i, q_j, w_i, w_j, LHS, RHS, vsign_i, vsign_j, xsign, qab_i, qab_j
        integer :: i, j
        real, parameter :: sigma = 2.0/3.0, alpha = 1.0, beta = 2.0

        a(1:n) = 0.0 
        do i=1, n
        do j=1, n+n_ghost
            if (i==j) cycle !skip itself (one less set of calculations)
            
            dist = abs(x(i)-x(j))
            xsign = sign(1.0, x(i)-x(j))

            q_i = dist/h(i)
            q_j = dist/h(j)

            call grad_cubic_spline(q_i, w_i)
            call grad_cubic_spline(q_j, w_j)

            select case (fix_choice)
            case (0,1) !no visocisty
            qab_i = 0
            qab_j = 0
            
            case (2) ! viscosity
            vsign_i = alpha*cs(i) - beta*(v(i)-v(j))*xsign !sign as in 1D so unit vector is scalar
            vsign_j = alpha*cs(j) - beta*(v(i)-v(j))*xsign !sign as in 1D so unit vector is scalar
                
            if ( (v(i)-v(j))*xsign < 0.0 ) then
                qab_i = -0.5 * rho(i)*vsign_i * (v(i)-v(j)) * xsign
                qab_j = -0.5 * rho(j)*vsign_j * (v(i)-v(j)) * xsign
            else 
                qab_i = 0
                qab_j = 0
            end if
            end select

            LHS = (pre(i) + qab_i) / rho(i)**2  * sigma/h(i)**2 * w_i * xsign
            RHS = (pre(j) + qab_j) / rho(j)**2  * sigma/h(j)**2 * w_j * xsign
                
            !Combining
            a(i) = a(i) - m(j) * (LHS + RHS)
                
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