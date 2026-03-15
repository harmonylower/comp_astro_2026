program SPH
    use hello, only:say
    use setup, only:init_iso_wave, init_shock_tube
    use write_output, only:output
    use calcs, only:get_den
    use boundary, only:set_ghosts
    use leapfrog, only:get_derivs, step
    use energy_calcs, only:write_ke_tot, get_ke, get_period

    implicit none
    integer, parameter :: nmax=250
    real, parameter :: dtout = 0.05, tmax = 5

    real, dimension(nmax):: x, v, m, h, rho, u, pre, cs, a
    real, dimension(int(tmax*1e4)) :: ke_tot, time_tot
    integer :: n, n_ghost, i
    real :: dt, tprint, time, xmax, xmin, dx

    integer :: problem = 1 !1:1D isothermal linear wave -- 2:isothermal shock tube -- 3:adiabatic eos
    integer :: fixes = 2!0:no fixes -- 1:+variable smoothing length -- 2:+artifical viscosity



    time = 0.0
    i=1

    !Inital calcuting and output
    call say()
    select case (problem)
    case (1) !isothermal linear wave
        xmax = 1.0
        xmin = 0.0
        call init_iso_wave(nmax, xmax, xmin, x, v, m, h, cs, n)
    case (2) !isothermal shock tube
        xmax = 0.5
        xmin = -0.5
        call init_shock_tube(x, rho, m, xmax, xmin, n)
    end select
    dx = xmax-xmin

    call set_ghosts(n, x, v, m, h, rho, u, pre, n_ghost, cs, dx, xmax,xmin)
    call get_derivs(n, n_ghost, x, m, h, rho, pre, cs, a, v, u, dt, problem,fixes, dx,xmax,xmin)
    call set_ghosts(n, x, v, m, h, rho, u, pre, n_ghost, cs, dx, xmax,xmin)

    print *, maxval(a)

    call output(time, x, v, a, m, h, rho, u, pre, n)

    tprint = dtout
    do while (time < tmax)
        call step(x, v, a, m, h, rho, pre, cs, dt, n, n_ghost, nmax, u, problem,fixes,dx,xmax,xmin)
        call set_ghosts(n, x, v, m, h, rho, u, pre, n_ghost,cs,dx,xmax,xmin)

        call get_ke(v(1:n), m(1:n), u(1:n), ke_tot)


        time = time + dt
        time_tot(i) = time
        i = i + 1

        if (time > tprint) then
            call output(time, x, v, a, m, h, rho, u, pre, n)
            tprint = tprint + dtout
        end if
    end do
    call write_ke_tot(ke_tot, time_tot)
    call get_period(ke_tot, time_tot)

end program SPH


