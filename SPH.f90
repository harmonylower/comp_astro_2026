program SPH
    use hello, only:say
    use setup, only:init
    use write_output, only:output
    use calcs, only:get_den
    use boundary, only:set_ghosts
    use leapfrog, only:get_derivs, step

    implicit none
    integer, parameter :: nmax=250
    real, dimension(nmax):: x, vel, mass, h, rho, u, pre, cs, acc
    integer :: n, n_ghost
    real :: time, dt, tprint

    real, parameter :: dtout = 0.05, tmax = 5

    time=0.0

    !Inital calcuting and output
    call say()
    call init(nmax, x, vel, mass, h, rho, pre, u, cs, n)
    call set_ghosts(n, x, vel, mass, h, rho, u, pre, n_ghost)
    call get_derivs(n, n_ghost, x, mass, h, rho, pre, cs, acc, vel, u)
    call output(time, x, vel, mass, h, rho, u, pre, n)

    dt = 0.5*h(1)/cs(1)
    print*, 'dt', dt
    tprint = dtout
    do while (time < tmax)
        call step(x, vel, acc, mass, h, rho, pre, cs, dt, n, n_ghost, nmax, u)
        call set_ghosts(n, x, vel, mass, h, rho, u, pre, n_ghost)
        time = time + dt
        if (time > tprint) then
            call output(time, x, vel, mass, h, rho, u, pre, n)
            tprint = tprint + dtout
        end if
    end do

end program SPH


