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
end module eos
