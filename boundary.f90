module boundary
contains
    subroutine set_ghosts(n, x, vel, mass, h, rho, u, pre, n_ghosts, cs, dx, xmax, xmin,problem)
        real, dimension(:), intent(inout) :: x, vel, mass, h, rho, u, pre, cs
        integer, intent(in) :: n, problem
        integer, intent(out) :: n_ghosts
        real, intent(in) :: dx, xmax, xmin
        integer :: i
        !PERIODIC BOUNDARY CONDITONS
        !PARTICLES FROM N+1 TO N+N_GHOST ARE GHOST PARTICLES
        
        n_ghosts = 0

        select case (problem)
        case (1) !isothermal linear sine wave
        do i=1,n
            !if 2 smothing lengths away is greater than the maximum xvalue 
            !then it must be include in the calculations for the particles with smaller x
            !so set a ghost particle with properties of the large x particle but a position before the low x particles
            if (x(i) + 2.0*h(i) > xmax) then  !the two is a parameter from the kernal, R kern = 2 for cubic spline 
                n_ghosts = n_ghosts + 1
                x(n+n_ghosts) = x(i) - dx !dx is total length not spacing
                vel(n+n_ghosts) = vel(i)
                mass(n+n_ghosts) = mass(i)
                h(n+n_ghosts) = h(i)
                rho(n+n_ghosts) = rho(i)
                u(n+n_ghosts) = u(i)
                pre(n+n_ghosts) = pre(i)
                cs(n+n_ghosts) = cs(i)

            elseif (x(i) - 2.0*h(i) < xmin) then
                n_ghosts = n_ghosts + 1
                x(n+n_ghosts) = x(i) + dx !dx is total length not spacing
                vel(n+n_ghosts) = vel(i)
                mass(n+n_ghosts) = mass(i)
                h(n+n_ghosts) = h(i)
                rho(n+n_ghosts) = rho(i)
                u(n+n_ghosts) = u(i)
                pre(n+n_ghosts) = pre(i)
                cs(n+n_ghosts) = cs(i)
            end if

        end do
    case (2) !fixed boundary conditions (shock tube and sod shock)
        nbound = 6
        vel(1:nbound) = 0.0
        vel(n-nbound:n) = 0.0
    end select
    end subroutine set_ghosts
end module boundary

