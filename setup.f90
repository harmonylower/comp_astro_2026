module init
    implicit none
    real, parameter :: xmax=1.0, xmin=0.0,  rho0=1.0
    
CONTAINS
    subroutine setup(nmax, x, vel, mass, h, rho, pre, u, cs, n)
        integer, intent(in) :: nmax
        integer, intent(out) :: n
        real, parameter :: PI=4.0*atan(1.0)
        real, dimension(nmax), intent(inout) :: x, vel, mass, h, rho, pre, cs, u
        real :: dx
        integer :: i

        dx = (xmax-xmin)/100.0 !ENSURE CHANGE LATER

        n = INT((xmax-xmin) / dx)

        do i=1,n
            cs(i)=1.0 !ENSURE CHANGED LATTER
            x(i) = xmin + real(i)*dx
            vel(i) = 1.0-4 * cs(i) * SIN((2.0*PI)/(xmax-xmin) * (x(i)-xmin))
            mass(i) = rho0 * dx 
            h(i) = 1.2 * dx
        end do



    end subroutine 
    
end module init

!put in safe guards for if n>nx