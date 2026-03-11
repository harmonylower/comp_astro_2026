module setup
    use eos, only: iso_eos
    real, parameter :: xmax=1.0, xmin=0.0,  rho0=1.0 !SOMEWHERE ELSE USES
    real, parameter :: hfac = 1.2 ! get_h uses    
    
contains
    subroutine init(nmax, x, vel, mass, h, rho, pre, u, cs, n, problem)
        integer, intent(in) :: nmax
        integer, intent(out) :: n
        real, parameter :: PI=4.0*atan(1.0)
        real, dimension(nmax), intent(inout) :: x, vel, mass, h, rho, pre, cs, u
        real :: dx
        integer :: i, problem

        select case (problem)
        case(1) !isothermal 1D linear wave

        dx = (xmax-xmin)/100.0 !ENSURE CHANGE LATER
        n = INT((xmax-xmin) / dx)
        do i=1,n
            cs(i)=1.0 !ENSURE CHANGED LATTER
            x(i) = xmin + (real(i)-0.5)*dx !starts particles with slight off set 
            vel(i) = 1.0E-4 * cs(i) * sin((2.0*PI)/(xmax-xmin) * (x(i)-xmin)) !ASK DANIEL
            mass(i) = rho0 * dx 
            h(i) = hfac * dx
            pre(i) = rho0 * cs(i)**2
        end do

        case(2) !isothermal shock tube problem

        end select
    end subroutine 
    
end module setup

!put in safe guards for if n>nx