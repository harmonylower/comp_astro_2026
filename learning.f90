program inital
    use hello, only:say
    IMPLICIT NONE
    INTEGER, PARAMETER :: nmax=3
    INTEGER :: i
    REAL, DIMENSION(nmax):: pos, vel, mass, sml, rho, u, pre, cs
    

    call say()


end program inital


!pi=4*atan(1.0)
!public is can use in any other place
!use say, only: say_hello