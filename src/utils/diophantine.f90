module diophantine
    use constants, only: dp

    implicit none
    private

    public :: gcd, lcm, rational_approx

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
        ! Best approximation x ~ p/q with q <= max_denom: the last continued
        ! fraction convergent or the largest semiconvergent that still fits
        real(dp), intent(in) :: x
        integer, intent(in) :: max_denom
        integer, intent(out) :: p, q
        integer :: a, k, p0, p1, q0, q1, pt, qt
        real(dp) :: r, xa

        p0 = 0; q0 = 1
        p1 = 1; q1 = 0

        if (max_denom <= 0) then
            error stop "max_denom must be positive"
        end if
        xa = abs(x)
        r = xa

        do
            a = int(r)
            pt = a*p1 + p0
            qt = a*q1 + q0
            if (qt > max_denom) then
                if (q1 + q0 <= max_denom) then
                    k = (max_denom - q0)/q1
                    pt = k*p1 + p0
                    qt = k*q1 + q0
                    if (abs(xa - real(pt, dp)/qt) < abs(xa - real(p1, dp)/q1)) then
                        p1 = pt
                        q1 = qt
                    end if
                end if
                exit
            end if
            p0 = p1; q0 = q1
            p1 = pt; q1 = qt
            if (abs(r - a) < 1.0d-12) exit
            r = 1.0_dp/(r - a)
        end do

        p = merge(p1, -p1, x >= 0.0_dp); q = q1
    end subroutine rational_approx

end module diophantine
