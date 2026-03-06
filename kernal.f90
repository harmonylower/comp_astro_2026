module kernal
    implicit none
contains
    subroutine cubic_spline(q,w)
        real, intent(in) :: q
        real, intent(out) :: w

        if (q >= 0.0 .and. q <1.0) then
            w = 0.25 * (2.0-q)**3 - (1.0-q)**3
        elseif (q >=1.0 .and. q<2.0) then
            w = 0.25 * (2.0-q)**3
        else
            w = 0.0
        end if
    end subroutine cubic_spline

    subroutine grad_cubic_spline(q, dw)
        real, intent(in) :: q
        real, intent(out) :: dw

        if (q >= 0.0 .and. q <1.0) then
            dw = -0.75 * (2.0-q)**2 + 3.0 * (1.0-q)**2
        elseif (q >=1.0 .and. q<2.0) then
            dw = -0.75 * (2.0-q)**2
        else
            dw = 0.0
        end if
    end subroutine grad_cubic_spline
end module kernal