module eos
    use setup, only:gamma
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
        integer :: i
  !for ideal gas gamma=5/3, for radiation pressure 4/2

        pre(1:n) = (gamma - 1) * rho(1:n) * u(1:n)
        
        do i = 1, n
            if (rho(i) > 0.0 .and. pre(i) > 0.0) then
                cs(i) = sqrt(gamma * pre(i) / rho(i))
            else
                cs(i) = 1.0  ! fallback value to prevent NaN
            end if
        end do

    end subroutine adiabatic_eos 
end module eos