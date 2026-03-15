module eos
    implicit none
contains
    subroutine iso_eos(n, rho, pre, cs)
        !isothermic equation of state
        real, dimension(:), intent(in) :: rho
        real, dimension(:), intent(inout) :: pre, cs
        integer, intent(in) :: n

        cs(1:n) = 1.0 

        pre(1:n) = rho(1:n) * cs(1:n)**2
    end subroutine iso_eos

    subroutine adiabatic_eos(n, rho, u, pre, cs)
        !Adiabatic equation of state
        real, dimension(:), intent(in) :: rho, u
        real, dimension(:), intent(inout) :: pre, cs
        integer, intent(in) :: n
        real, parameter :: gamma = 5.0/3.0 !for ideal gas, for radiation pressure 4/2


        pre(1:n) = (gamma - 1) * rho(1:n) * u(1:n)
        cs(1:n) = sqrt(gamma * Pre/rho)
    end subroutine adiabatic_eos 
end module eos