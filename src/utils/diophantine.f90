module diophantine
    use constants, only: dp

    implicit none

contains

    function gcd(a, b) result(g)
        integer, intent(in) :: a, b
        integer :: g, aa, bb, tmp

        aa = abs(a); bb = abs(b)
        do while (bb /= 0)
            tmp = bb; bb = mod(aa, bb); aa = tmp
        end do
        g = aa
    end function gcd

    function lcm(a, b) result(l)
        integer, intent(in) :: a, b
        integer :: l

        l = abs(a/gcd(a, b)*b)   ! divide first to avoid overflow
    end function lcm

    subroutine rational_approx(x, max_denom, p, q)
        ! Approximate x ~ p/q via continued fraction convergents, q <= max_denom
        real(dp), intent(in) :: x
        integer, intent(in) :: max_denom
        integer, intent(out) :: p, q
        integer :: a, p0, p1, q0, q1, pt, qt
        real(dp) :: r

        p0 = 0; q0 = 1
        p1 = 1; q1 = 0

        if (max_denom <= 0) then
            error stop "max_denom must be positive"
        end if
        r = abs(x)

        do
            a = int(r)
            pt = a*p1 + p0
            qt = a*q1 + q0
            if (qt > max_denom) exit
            p0 = p1; q0 = q1
            p1 = pt; q1 = qt
            if (abs(r - a) < 1.0d-12) exit
            r = 1.0_dp/(r - a)
        end do

        p = merge(p1, -p1, x >= 0.0_dp); q = q1
    end subroutine rational_approx

end module diophantine
