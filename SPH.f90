program SPH
    use hello, only:say
    use init, only:setup
    use write_output, only:output
    use get_density, only:get_den
    use boundary, only:set_ghosts

    implicit none
    integer, parameter :: nmax=140
    real, dimension(nmax):: x, vel, mass, h, rho, u, pre, cs
    integer :: n, n_ghost
    real :: time

    time=0.0
    call say()
    call setup(nmax, x, vel, mass, h, rho, pre, u, cs, n)
    
    call set_ghosts(n, x, vel, mass, h, rho, u, pre, n_ghost)

    call get_den(n, n_ghost, x, mass, h, rho)


    call output(time, x, vel, mass, h, rho, u, pre, n)
end program SPH

!parameter from the kernal, R kern = 2 for cubic spline hence 2*smoothing length is criteria
!split kernal and density into two seperate so that v and u can be used in same module as density
!use functinos
!get rid of kind 8 and d0 on numbers
!consider where sound speed is coming from when putting in module
!no shouting :(
!kernal out