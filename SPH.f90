program SPH
    use hello, only:say
    use setup, only:init, hfac
    use write_output, only:output
    use calcs, only:get_den
    use boundary, only:set_ghosts
    use leapfrog, only:get_derivs, step
    use energy_calcs, only:write_ke_tot, get_ke, get_period

    implicit none
    integer, parameter :: nmax=250
    real, parameter :: dtout = 0.05, tmax = 5

    real, dimension(nmax):: x, v, m, h, rho, u, pre, cs, a
    real, dimension(int(tmax*1e4)) :: ke_tot
    integer :: n, n_ghost, i
    real :: time, dt, tprint

    integer :: problem = 1 !1:1D isothermal linear wave -- 2:isothermal shock tube 
    integer :: fixes = 0 !0:no fixes -- 1:+variable smoothing length -- 2:+artifical viscosity



    time=0.0

    !Inital calcuting and output
    call say()
    call init(nmax, x, v, m, h, rho, pre, u, cs, n,problem)
    call set_ghosts(n, x, v, m, h, rho, u, pre, n_ghost,cs)
    call get_derivs(n, n_ghost, x, m, h, rho, pre, cs, a, v, u)

    call output(time, x, v, a, m, h, rho, u, pre, n)
    
    dt = 0.25*minval(h(1:n)/cs(1:n)) !ENSURE CHANGE LATER

    tprint = dtout
    do while (time < tmax)
        call step(x, v, a, m, h, rho, pre, cs, dt, n, n_ghost, nmax, u)
        call set_ghosts(n, x, v, m, h, rho, u, pre, n_ghost,cs)

        call get_ke(v(1:n), m(1:n), u(1:n), ke_tot)


        time = time + dt
        if (time > tprint) then
            call output(time, x, v, a, m, h, rho, u, pre, n)
            tprint = tprint + dtout
        end if
    end do
    call write_ke_tot(ke_tot, tmax, dt)
    call get_period(ke_tot, dt)

end program SPH


