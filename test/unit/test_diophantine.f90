program test_diophantine
    implicit none

    logical :: test_failed

    test_failed = .false.

    call test_gcd(test_failed)
    call test_lcm(test_failed)
    call test_rational_approx(test_failed)

    if (test_failed) error stop

contains

    subroutine test_gcd(test_failed)
        use diophantine, only: gcd

        logical, intent(inout) :: test_failed

        if (gcd(12, 8) /= 4) then
            print *, "-------------------------------------------------------------"
            print *, "test_diophantine failed: gcd(12, 8) /= 4"
            test_failed = .true.
        end if

        if (gcd(7, 13) /= 1) then
            print *, "-------------------------------------------------------------"
            print *, "test_diophantine failed: gcd(7, 13) /= 1"
            test_failed = .true.
        end if

        if (gcd(0, 5) /= 5) then
            print *, "-------------------------------------------------------------"
            print *, "test_diophantine failed: gcd(0, 5) /= 5"
            test_failed = .true.
        end if

        if (gcd(5, 0) /= 5) then
            print *, "-------------------------------------------------------------"
            print *, "test_diophantine failed: gcd(5, 0) /= 5"
            test_failed = .true.
        end if

        if (gcd(100, 75) /= 25) then
            print *, "-------------------------------------------------------------"
            print *, "test_diophantine failed: gcd(100, 75) /= 25"
            test_failed = .true.
        end if

        if (gcd(-12, 8) /= 4) then
            print *, "-------------------------------------------------------------"
            print *, "test_diophantine failed: gcd(-12, 8) /= 4"
            test_failed = .true.
        end if
    end subroutine test_gcd

    subroutine test_lcm(test_failed)
        use diophantine, only: lcm

        logical, intent(inout) :: test_failed

        if (lcm(4, 6) /= 12) then
            print *, "-------------------------------------------------------------"
            print *, "test_diophantine failed: lcm(4, 6) /= 12"
            test_failed = .true.
        end if

        if (lcm(3, 7) /= 21) then
            print *, "-------------------------------------------------------------"
            print *, "test_diophantine failed: lcm(3, 7) /= 21"
            test_failed = .true.
        end if

        if (lcm(12, 18) /= 36) then
            print *, "-------------------------------------------------------------"
            print *, "test_diophantine failed: lcm(12, 18) /= 36"
            test_failed = .true.
        end if

        if (lcm(5, 5) /= 5) then
            print *, "-------------------------------------------------------------"
            print *, "test_diophantine failed: lcm(5, 5) /= 5"
            test_failed = .true.
        end if
    end subroutine test_lcm

    subroutine test_rational_approx(test_failed)
        use constants, only: dp, pi
        use diophantine, only: rational_approx

        logical, intent(inout) :: test_failed
        integer :: p, q

        call rational_approx(-pi, 1000, p, q)
        if (p /= -355 .or. q /= 113) then
            print *, "-------------------------------------------------------------"
            print *, "test_diophantine failed: rational_approx(pi, 1000)"
            print *, "found: ", p, "/", q
            print *, "expected: -355 / 113"
            test_failed = .true.
        end if

        call rational_approx(1.0_dp/3.0_dp, 100, p, q)
        if (p /= 1 .or. q /= 3) then
            print *, "-------------------------------------------------------------"
            print *, "test_diophantine failed: rational_approx(1/3, 100)"
            print *, "found: ", p, "/", q
            print *, "expected: 1 / 3"
            test_failed = .true.
        end if

        call rational_approx(0.75_dp, 100, p, q)
        if (p /= 3 .or. q /= 4) then
            print *, "-------------------------------------------------------------"
            print *, "test_diophantine failed: rational_approx(0.75, 100)"
            print *, "found: ", p, "/", q
            print *, "expected: 3 / 4"
            test_failed = .true.
        end if

        ! semiconvergent 140/99 is closer than the convergent 99/70
        call rational_approx(sqrt(2.0_dp), 100, p, q)
        if (p /= 140 .or. q /= 99) then
            print *, "-------------------------------------------------------------"
            print *, "test_diophantine failed: rational_approx(sqrt(2), 100)"
            print *, "found: ", p, "/", q
            print *, "expected: 140 / 99"
            test_failed = .true.
        end if

        ! the next convergent 1/9 does not fit; 1/8 is closer than 0/1
        call rational_approx(0.1094_dp, 8, p, q)
        if (p /= 1 .or. q /= 8) then
            print *, "-------------------------------------------------------------"
            print *, "test_diophantine failed: rational_approx(0.1094, 8)"
            print *, "found: ", p, "/", q
            print *, "expected: 1 / 8"
            test_failed = .true.
        end if
    end subroutine test_rational_approx

end program test_diophantine
