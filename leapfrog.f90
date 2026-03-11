module leapfrog
    use calcs, only: get_accel, get_den, get_h
    use eos, only: iso_eos
    use boundary, only:set_ghosts
    implicit none
    
contains
    subroutine get_derivs(n, n_ghost, x, m, h, rho, pre, cs, a, v,u)
        integer, intent(in) :: n
        integer, intent(inout) :: n_ghost
        real, dimension(:), intent(inout) :: x, m, h, u
        real, dimension(:), intent(inout) :: rho, pre, cs, a, v
        integer :: i

        do i=1, 3
            call set_ghosts(n, x, v, m, h, rho, u, pre, n_ghost,cs) 
            call get_den(n, n_ghost, x, m, h, rho)
            call get_h(m, rho, n, h)
        end do
        call set_ghosts(n, x, v, m, h, rho, u, pre, n_ghost,cs) !remove once work out what initial pressure should be
        call iso_eos(n, rho, pre, cs)
        call set_ghosts(n, x, v, m, h, rho, u, pre, n_ghost,cs) !remove once work out what initial pressure should be
        call get_accel(n, n_ghost, x, v, m, h, rho, pre, a, cs)
    
    
    end subroutine get_derivs
    
    subroutine step(x, v, a, m, h, rho, pre, cs, dt, n, n_ghost,nmax,u)
        real, dimension(:), intent(inout) :: x, v, a, m, h, rho, pre, cs, u
        real, intent(in) :: dt
        integer, intent(in) :: n, nmax
        integer, intent(inout) :: n_ghost
        real, dimension(nmax) :: a0, v_half
        a0 = a
        x = x + dt*v + 0.5 * dt**2 * a0
        v_half = v + dt*a0
        call get_derivs(n, n_ghost, x, m, h, rho, pre, cs, a, v_half, u)
        v = v_half + 0.5*dt* (a - a0)
    end subroutine step
    
end module leapfrog