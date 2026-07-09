program test_integrate_1d_substituted
    use constants, only: dp, pi
    use integrate_substituted, only: integrate_1d_substituted
    use utils, only: not_same

    implicit none

    real(dp), parameter :: reltol = 1.0e-8_dp, abstol = 0.0_dp
    real(dp), dimension(2), parameter :: interval_0 = [-1.0_dp, 1.0_dp]
    real(dp), parameter :: integral_0 = 2.0_dp/3.0_dp
    real(dp), dimension(2), parameter :: interval_1 = [-1.0_dp, 1.0_dp]
    real(dp), parameter :: integral_1 = 0.5_dp*pi
    real(dp), dimension(2), parameter :: interval_2 = [0.0_dp, 1.0_dp]
    real(dp), parameter :: integral_2 = 2.0_dp/3.0_dp
    real(dp), dimension(2), parameter :: interval_3 = [0.0_dp, 0.5_dp*pi]
    real(dp), parameter :: integral_3 = 2.0_dp*(sqrt(2.0_dp) - 1.0_dp)
    real(dp), dimension(2), parameter :: interval_4 = [-pi, pi]
    real(dp), parameter :: integral_4 = 4.0_dp*sqrt(2.0_dp)
    real(dp), dimension(2), parameter :: interval_5 = [-pi, pi]
    real(dp), parameter :: eps = 0.5_dp
    real(dp), parameter :: integral_5 = integral_4*(1.0_dp + eps/3.0_dp)

    real(dp) :: found_integral
    logical :: failed_test

    failed_test = .false.

    call integrate_1d_substituted(trial_func_0, &
                                  interval_0(1), &
                                  interval_0(2), &
                                  found_integral)
    if (not_same(integral_0, found_integral, reltol, abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_integrate_1d_substituted failed: integral 0"
        print *, "found: ", found_integral
        print *, "analytic: ", integral_0
        failed_test = .true.
    end if

    call integrate_1d_substituted(trial_func_1, &
                                  interval_1(1), &
                                  interval_1(2), &
                                  found_integral)
    if (not_same(integral_1, found_integral, reltol, abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_integrate_1d_substituted failed: integral 1"
        print *, "found: ", found_integral
        print *, "analytic: ", integral_1
        failed_test = .true.
    end if

    call integrate_1d_substituted(trial_func_2, &
                                  interval_2(1), &
                                  interval_2(2), &
                                  found_integral)
    if (not_same(integral_2, found_integral, reltol, abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_integrate_1d_substituted failed: integral 2"
        print *, "found: ", found_integral
        print *, "analytic: ", integral_2
        failed_test = .true.
    end if

    call integrate_1d_substituted(trial_func_3, &
                                  interval_3(1), &
                                  interval_3(2), &
                                  found_integral)
    if (not_same(integral_3, found_integral, reltol, abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_integrate_1d_substituted failed: integral 3"
        print *, "found: ", found_integral
        print *, "analytic: ", integral_3
        failed_test = .true.
    end if

    call integrate_1d_substituted(trial_func_4, &
                                  interval_4(1), &
                                  interval_4(2), &
                                  found_integral)
    if (not_same(integral_4, found_integral, reltol, abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_integrate_1d_substituted failed: integral 4"
        print *, "found: ", found_integral
        print *, "analytic: ", integral_4
        print *, "relative error: ", found_integral/integral_4 - 1.0_dp
        failed_test = .true.
    end if

    call integrate_1d_substituted(trial_func_5, &
                                  interval_5(1), &
                                  interval_5(2), &
                                  found_integral)
    if (not_same(integral_5, found_integral, reltol, abstol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_integrate_1d_substituted failed: integral 5"
        print *, "found: ", found_integral
        print *, "analytic: ", integral_5
        print *, "relative error: ", found_integral/integral_5 - 1.0_dp
        failed_test = .true.
    end if

    if (failed_test) error stop

contains
    function trial_func_0(x)
        real(dp), intent(in) :: x
        real(dp) :: trial_func_0

        trial_func_0 = x**2
    end function trial_func_0

    function trial_func_1(x)
        real(dp), intent(in) :: x
        real(dp) :: trial_func_1

        trial_func_1 = sqrt(1.0_dp - x**2)
    end function trial_func_1

    function trial_func_2(x)
        real(dp), intent(in) :: x
        real(dp) :: trial_func_2

        trial_func_2 = sqrt(1.0_dp - x)
    end function trial_func_2

    function trial_func_3(x)
        real(dp), intent(in) :: x
        real(dp) :: trial_func_3

        trial_func_3 = sqrt(1.0_dp - sin(x))
    end function trial_func_3

    function trial_func_4(x)
        real(dp), intent(in) :: x
        real(dp) :: trial_func_4

        trial_func_4 = sqrt(1.0_dp + cos(x))
    end function trial_func_4

    function trial_func_5(x)
        real(dp), intent(in) :: x
        real(dp) :: trial_func_5

        trial_func_5 = trial_func_4(x)*(1.0_dp + eps*cos(x))
    end function trial_func_5

end program test_integrate_1d_substituted
