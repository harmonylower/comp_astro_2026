module leapfrog
    use calcs, only: get_accel, get_den, get_h
    use eos, only: iso_eos, adiabatic_eos
    use boundary, only:set_ghosts
    use setup, only: fix_choice, init_choice, eos_choice
    implicit none
    
contains
    subroutine get_derivs(n, n_ghost, x, m, h, rho, pre, cs, a, v,u,du, dt_new, dx, xmax,xmin)
        integer, intent(in) :: n
        integer, intent(inout) :: n_ghost
        real, dimension(:), intent(inout) :: x, m, h, u, du
        real, dimension(:), intent(inout) :: rho, pre, cs, a, v
        real, intent(out) :: dt_new
        real, intent(in) :: dx, xmax, xmin
        integer :: i
        
        select case (fix_choice)
        case (0)
            call set_ghosts(n, x, v, m, h, rho, u, pre, n_ghost, cs, dx, xmax, xmin) 
            call get_den(n, n_ghost, x, m, h, rho)
        case (1,2)
        do i=1, 3
            call set_ghosts(n, x, v, m, h, rho, u, pre, n_ghost,cs, dx, xmax,xmin) 
            call get_den(n, n_ghost, x, m, h, rho)
            call get_h(m, rho, n, h)
        end do
        end select
        
        call set_ghosts(n, x, v, m, h, rho, u, pre, n_ghost,cs, dx, xmax,xmin) !remove once work out what initial pressure should be
        
        select case (eos_choice)
        case (1) !isothermal eos
            call iso_eos(n, rho, pre, cs)
        case (2) !adiabatic eos
            call adiabatic_eos(n, rho, u, pre, cs)
        end select

    
        call set_ghosts(n, x, v, m, h, rho, u, pre, n_ghost,cs, dx, xmax, xmin) !remove once work out what initial pressure should be
        call get_accel(n, n_ghost, x, v, m, h, rho, pre, a, cs, du)
    
        dt_new = 0.25*minval(h(1:n)/cs(1:n)) 

    
    end subroutine get_derivs
    
    subroutine step(x, v, a, m, h, rho, pre, cs, dt, n, n_ghost,nmax,u,du,dx,xmax,xmin)
        real, dimension(:), intent(inout) :: x, v, a, m, h, rho, pre, cs, u, du
        real, intent(inout) :: dt
        real, intent(in) :: dx, xmax, xmin
        integer, intent(in) :: n, nmax
        integer, intent(inout) :: n_ghost
        real, dimension(nmax) :: a0, v_half, du0, u_half
        real :: dt_new
        a0 = a
        du0 = du
        x = x + dt*v + 0.5 * dt**2 * a0
        v_half = v + dt*a0
        u_half = u + dt*du0
        call get_derivs(n, n_ghost, x, m, h, rho, pre, cs, a, v_half, u, du, dt_new, dx, xmax,xmin)
        v = v_half + 0.5*dt* (a - a0)
        u = u_half + 0.5*dt* (du - du0)
        dt = dt_new
    end subroutine step
    
end module leapfrog