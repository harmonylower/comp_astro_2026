module write_output
    implicit none
    integer :: nfile = -1
    
contains
    subroutine output(time, x, vel, mass, h, rho, u, pre, n)
        implicit none
        integer, intent(in) :: n
        real, dimension(:), intent(in) :: x, vel, mass, h, rho, pre, u
        real, intent(in) :: time
        
        integer :: i, lu
        character(len=256) :: filename

        nfile = nfile + 1
        !WRITES PARTICLE INFORMATION TO A FILE. EXCLUDING GHOST PARTICLES
        write(filename,"(a,i5.5,a)") 'output_',nfile,'.txt'
        open(newunit=lu, file=filename, status='replace', action='write')
        
        write(lu,*) '# x, vel, mass, smoothing length, density, internal energy, pressure, sound speed'
        write(lu,*) time
        do i=1,n
            write(lu,*) x(i), vel(i), mass(i), h(i), rho(i), u(i), pre(i)
        enddo

        close(lu)
    end subroutine output



end module write_output


!rewrite to take parameters when naming