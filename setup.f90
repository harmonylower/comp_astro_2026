module setup
    real, parameter :: hfac = 1.2, gamma=1.4 ! get_h uses 
    integer :: init_choice = 3 !1:1D isothermal linear wave -- 2:shock tube -- 3:sod shock (also sets boundary condition)
    integer :: eos_choice = 2 !1:isothermic 2:adiabaitc
    integer :: fix_choice = 2!0:no fixes -- 1:+variable smoothing length -- 2:+artifical viscosity
   
    
contains
    subroutine init_iso_wave(nmax, xmax, xmin, x, vel, mass, h, cs, n)
        integer, intent(in) :: nmax
        integer, intent(out) :: n
        real, parameter :: PI=4.0*atan(1.0), rho0=1.0
        real, dimension(nmax), intent(inout) :: x, vel, mass, h, cs
        integer :: i
        real :: x_step

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
    
    subroutine init_shock_tube(x, v, rho, mass, h, xmax, xmin, n)
        integer, intent(out) :: n
        real, dimension(:), intent(inout) :: rho, x, v, mass, h
        real, intent(in) :: xmax, xmin
        real :: x_step, xval
        n = 0
        xval = xmin
        do while (xval<=xmax)
            n = n+1
            x(n) = xval

            if (xval<0) then
                rho(n) = 1.0
                x_step = 0.001
                mass(n) = x_step * rho(n)
                h(n) = hfac * x_step

            else
                rho(n) = 0.1
                x_step = 0.01
                mass(n) = x_step * rho(n)
                h(n) = hfac * x_step
            end if

            v(n)=0.0
            xval = xval + x_step
        end do

    end subroutine init_shock_tube

subroutine init_sod_shock(x, v, rho, pre, mass, h, u, xmax, xmin, n)
        integer, intent(out) :: n
        real, dimension(:), intent(inout) :: rho, x, v, mass, h, pre, u
        real, intent(in) :: xmax, xmin
        real :: x_step, xval
        n = 0
        xval = xmin
        do while (xval<=xmax)
            n = n+1
            x(n) = xval

            if (xval<0) then
                rho(n) = 1.0
                pre(n) = 1.0
                x_step = 0.001
            else
                rho(n) = 0.125
                pre(n) = 0.1
                x_step = 0.008
            end if

            mass(n) = x_step * rho(n)
            h(n) = hfac * x_step
            v(n) = 0.0
            u(n) = pre(n) / (rho(n)*(gamma-1))
            xval = xval + x_step
        end do

    end subroutine init_sod_shock
end module setup

!put in safe guards for if n>nx
! need to output/initallise cs
!work out whats happening with velocity