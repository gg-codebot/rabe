program test_find_extrema
    use constants, only: dp, pi
    use utils, only: not_same

    implicit none

    call test_find_local_minima()
    call test_find_local_maxima()
    call test_find_global_extrema()
    call test_find_global_minimum()

contains

    subroutine test_find_local_minima()
        use find_extrema, only: find_local_minima

        ! as expect_locs is give with 5 significant digits
        real(dp), parameter :: tol = 1e-4
        real(dp), parameter :: error_limit = 1e-7
        real(dp) :: tol_2
        real(dp), dimension(2), parameter :: interval = [0.0_dp, 2.0_dp*pi]
        real(dp), dimension(2), parameter :: expect_locs = [2.1386_dp, 5.64832_dp]
        real(dp), dimension(2), parameter :: interval_2 = [-1.5_dp, 3.5_dp]
        real(dp), dimension(2), parameter :: expect_locs_2 = [-1.0_dp, 3.0_dp]
        real(dp), dimension(:), allocatable :: found_locs
        real(dp), dimension(:), allocatable :: achieved_error

        call find_local_minima(negative_sincos_func, interval, found_locs)

        if (size(found_locs) /= 2) then
            print *, "-------------------------------------------------------------"
            print *, "test_find_local_minima"
            print *, "found ", size(found_locs), " minima, expected 2"
            error stop
        end if

        if (not_same(expect_locs, found_locs, reltol_in=tol, abstol_in=0.0_dp)) then
            print *, "-------------------------------------------------------------"
            print *, "test_find_local_minima"
            print *, "found locations: ", found_locs
            print *, "expected locations: ", expect_locs
            error stop
        end if

        call find_local_minima(poly_func, interval_2, found_locs, achieved_error)
        if (maxval(achieved_error) > error_limit) then
            print *, "-------------------------------------------------------------"
            print *, "test_find_local_minima"
            print *, "achieved error is too large: ", achieved_error
            print *, "error limit: ", error_limit
            error stop
        end if

        tol_2 = maxval(achieved_error)*2.0_dp
        if (not_same(expect_locs_2, found_locs, abstol_in=tol_2)) then
            print *, "-------------------------------------------------------------"
            print *, "test_find_local_minima"
            print *, "found locations: ", found_locs
            print *, "expected locations: ", expect_locs_2
            print *, "actual error: ", abs(found_locs - expect_locs_2)
            print *, "diagnosed achieved error: ", achieved_error
            error stop
        end if
    end subroutine test_find_local_minima

    subroutine negative_sincos_func(x, value)
        real(dp), dimension(:), intent(in) :: x
        real(dp), dimension(:), intent(out) :: value

        call sincos_func(x, value)
        value = -value
    end subroutine negative_sincos_func

    subroutine poly_func(x, value)
        real(dp), dimension(:), intent(in) :: x
        real(dp), dimension(:), intent(out) :: value

        value = 0.25_dp*x**4.0_dp - 2.0_dp/3.0_dp*x**3.0_dp - 1.5_dp*x**2.0_dp
    end subroutine poly_func

    subroutine test_find_local_maxima()
        use find_extrema, only: find_local_maxima

        real(dp), parameter :: tol = 1e-4
        real(dp), dimension(2), parameter :: interval = [0.0_dp, 2.0_dp*pi]
        real(dp), dimension(2), parameter :: expect_locs = [2.1386_dp, 5.64832_dp]
        real(dp), dimension(:), allocatable :: found_locs

        call find_local_maxima(sincos_func, interval, found_locs)

        if (size(found_locs) /= 2) then
            print *, "-------------------------------------------------------------"
            print *, "test_find_local_maxima"
            print *, "found ", size(found_locs), " maxima, expected 2"
            error stop
        end if

        if (any(abs(found_locs/expect_locs - 1.0_dp) > tol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_find_local_maxima"
            print *, "found locations: ", found_locs
            print *, "expected locations: ", expect_locs
            error stop
        end if

    end subroutine test_find_local_maxima

    subroutine sincos_func(x, value)
        real(dp), dimension(:), intent(in) :: x
        real(dp), dimension(:), intent(out) :: value

        value = cos(x) - sin(2*x)
    end subroutine sincos_func

    subroutine test_find_global_extrema
        use find_extrema, only: find_global_extrema

        real(dp) :: extrema(2), expected_extrema(2)
        real(dp), parameter :: interval(2) = [0.0_dp, 2.0_dp*pi]
        real(dp), parameter :: range = abs(interval(2) - interval(1))
        real(dp), parameter :: tol = 1e-3

        expected_extrema = [0.5_dp, 1.5_dp]*pi
        extrema = find_global_extrema(sin_func, interval)

        if (any(abs(extrema - expected_extrema)/range > tol)) then
            print *, "-------------------------------------------------------------"
            print *, "test_find_global_extrema"
            print *, "found extrema: ", extrema
            print *, "expected extrema: ", expected_extrema
            print *, "relative error: ", abs(extrema - expected_extrema)/range
            error stop
        end if
    end subroutine test_find_global_extrema

    subroutine test_find_global_minimum
        use find_extrema, only: find_global_minimum

        real(dp) :: min_loc, expected
        real(dp), parameter :: interval(2) = [0.0_dp, 2.0_dp*pi]
        real(dp), parameter :: reltol = 1e-3

        expected = 1.5_dp*pi
        min_loc = find_global_minimum(sin_func, interval, reltol)

        if (not_same(min_loc, expected, reltol_in=reltol, abstol_in=0.0_dp)) then
            print *, "-------------------------------------------------------------"
            print *, "test_find_global_minimum"
            print *, "found location: ", min_loc
            print *, "expected location: ", expected
            error stop
        end if
    end subroutine test_find_global_minimum

    subroutine sin_func(x, value)
        real(dp), dimension(:), intent(in) :: x
        real(dp), dimension(:), intent(out) :: value

        value = sin(x)
    end subroutine sin_func

end program test_find_extrema
