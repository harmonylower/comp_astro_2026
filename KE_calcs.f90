module energy_calcs
    implicit none
    integer :: entryno=1
    contains
    subroutine get_ke(v, m, ke, ke_tot)
        real, dimension(:), intent(in) :: v, m
        real, dimension (:), intent(inout) :: ke, ke_tot
        ke = 0.5 * m * v**2
       
        ke_tot(entryno) = sum(ke)
        entryno = entryno + 1
    end subroutine get_ke

    subroutine write_ke_tot(ke_tot, tmax, dt)
        real, dimension(:), intent(in) :: ke_tot
        real, intent(in) :: tmax, dt
        integer :: i, lu

        open(lu, file='ke_tot.txt', status='replace', action='write')
        write(lu,*) '# time (s), Total Kinetic Energy'
        
        do i = 1, int(tmax/dt)
            write(lu,*) i*dt, ke_tot(i)
        end do
        close(lu)
    end subroutine write_ke_tot

    subroutine get_period(ke_tot, dt)
        real, dimension(:), intent(in) :: ke_tot
        real, intent(in) :: dt
        real :: period
        integer :: i, n

        n = size(ke_tot)
        do i=2,n-1
            if (ke_tot(i-1) < ke_tot(i) .and. ke_tot(i+1) < ke_tot(i)) then
                period = 2.0 * i*dt !2*distance between 2 peaks
                print *, 'Period = ', period, 's'
                return !only finds distance for first peak
            end if
        end do

    end subroutine get_period


end module energy_calcs