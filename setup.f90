module setup
    use eos, only: iso_eos
    !SOMEWHERE ELSE USES xmax and xmin
    real, parameter :: hfac = 1.2 ! get_h uses    
    
contains
    subroutine init_iso_wave(nmax, xmax, xmin, x, vel, mass, h, cs, n)
        integer, intent(in) :: nmax
        integer, intent(out) :: n
        real, parameter :: PI=4.0*atan(1.0), rho0=1.0
        real, dimension(nmax), intent(inout) :: x, vel, mass, h, cs
        integer :: i

        x_step = (xmax-xmin)/100.0
        n = int((xmax-xmin) / x_step)
        do i=1,n
            cs(i)=1.0 !ENSURE CHANGED LATTER
            x(i) = xmin + (real(i)-0.5)*x_step !starts particles with slight off set 
            vel(i) = 1.0E-4 * cs(i) * sin((2.0*PI)/(xmax-xmin) * (x(i)-xmin)) !ASK DANIEL
            mass(i) = rho0 * x_step 
            h(i) = hfac * x_step
        end do
    end subroutine 
    
    subroutine init_shock_tube(x, rho, mass, xmax, xmin, n)
        integer, intent(out) :: n
        real, dimension(:), intent(inout) :: rho, x, mass
        real, intent(in) :: xmax, xmin
        real :: x_step
        n = 1
        xval = xmin
        do while (xval<=xmax)
            x(n) = xval
            if (xval<0) then
                rho(n) = 1.0
                x_step = 0.001
                mass(n) = x_step * rho(n)
            else
                rho(n) = 0.1
                x_step = 0.01
                mass(n) = x_step * rho(n)
            end if
            xval = xval + x_step
            n = n+1
        end do

    end subroutine init_shock_tube
end module setup

!put in safe guards for if n>nx