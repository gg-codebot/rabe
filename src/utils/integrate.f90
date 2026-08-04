module integrate
    use constants, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan

    implicit none

    integer, parameter :: quadkind = 8

    abstract interface
        function integrand_i(x)
            use constants, only: dp
            real(dp), intent(in) :: x
            real(dp) :: integrand_i
        end function integrand_i
    end interface

    procedure(integrand_i), private, pointer, save :: integrand => null()

contains

    subroutine integrate_1d(f, a, b, result)
        use quadpack_double, only: qag => dqag
        procedure(integrand_i) :: f
        real(dp), intent(in) :: a, b
        real(dp), intent(out) :: result

        real(quadkind), parameter :: abs_error_tol_quadkind = 0.0_dp
        real(quadkind), parameter :: rel_error_tol_quadkind = 1.0e-6
        integer, parameter :: order_key = 6
        integer, parameter :: limit = 500
        integer, parameter :: lenw = limit*4
        real(quadkind) :: a_quadkind, b_quadkind, result_quadkind
        real(quadkind) :: abs_error, error_limit
        integer :: error_flag

        logical :: error_occurred

        integer :: neval_dummy, last
        integer :: iwork(limit)
        real(quadkind) :: work(lenw)

        if (ieee_is_nan(a) .or. ieee_is_nan(b)) then
            print *, "Integration limits must not be NaN!"
            print *, "a =", a, "b =", b
            error stop
        end if

        ! quadpack operates in and requires real(8)
        a_quadkind = convert_to_quadkind(a)
        b_quadkind = convert_to_quadkind(b)

        integrand => f
        ! integration by globally adaptive interval subdivision (quadpack import)
        call qag(quadkind_integrand, &
                 a_quadkind, &
                 b_quadkind, &
                 abs_error_tol_quadkind, &
                 rel_error_tol_quadkind, &
                 order_key, &
                 result_quadkind, &
                 abs_error, &
                 neval_dummy, &
                 error_flag, &
                 limit, &
                 lenw, &
                 last, &
                 iwork, &
                 work)
        integrand => null()

        error_limit = abs(result_quadkind)*rel_error_tol_quadkind &
                      + abs_error_tol_quadkind
        error_occurred = .false.

        if (abs_error > error_limit) then
            print *, "Integration warning: absolute error =", abs_error
            print *, "bigger than required ", error_limit
            print *, "relative error", abs_error/abs(result_quadkind)
            error_occurred = .true.
        end if

        if (error_flag /= 0) then
            print *, "Integration warning: error =", error_flag
            error_occurred = .true.
        end if

        result = convert_to_dp(result_quadkind)

        if (ieee_is_nan(result)) then
            print *, "Integration result is NaN!"
            error_occurred = .true.
        end if

        if (error_occurred) error stop

    end subroutine integrate_1d

    !> The integrand for quadpack needs to be of the below signature.
    !> Note for developer: due to quadpack interface, functions with parameters
    !> $$
    !> f(x, p_1, p_2, ...)
    !> $$
    !> need to be wrapped
    !>$$
    !>f(x, p_1, p_2, ...) \to f(x)
    !>$$
    !> which at some point in the call chain requires access
    !> to outer scopes where parameters are defined. If one wants to avoid
    !> executable stacks, one has to resort to module variables for that.
    !> Type bound procedures which carry the parameters with them fail latest
    !> at this boundary here. The call chain can not be made thread safe.
    !> We opt therefore for the most straight forward solution i.e. global state.
    function quadkind_integrand(x_quadkind)
        real(quadkind), intent(in) :: x_quadkind
        real(quadkind) :: quadkind_integrand

        real(dp) :: x_dp

        x_dp = real(x_quadkind, kind=dp)
        quadkind_integrand = real(integrand(x_dp), kind=quadkind)
    end function quadkind_integrand

    function convert_to_quadkind(val_dp) result(val_quadkind)
        real(dp), intent(in) :: val_dp
        real(quadkind) :: val_quadkind

        val_quadkind = real(val_dp, kind=quadkind)
    end function convert_to_quadkind

    function convert_to_dp(val_quadkind) result(val_dp)
        real(quadkind), intent(in) :: val_quadkind
        real(dp) :: val_dp

        val_dp = real(val_quadkind, kind=dp)
    end function convert_to_dp

    function sum_trapez_1d(x, y)
        real(dp), dimension(:), intent(in) :: x, y
        real(dp) :: sum_trapez_1d

        integer :: n

        n = size(x)
        if (.not. n == size(y)) then
            print *, "Error in sum_trapez_1d: x and y need to have same length!"
            error stop
        end if
        if (n < 2) then
            print *, "Error in sum_trapez_1d: len(x) must be at least 2!"
            error stop
        end if

        sum_trapez_1d = 0.5_dp*sum((x(2:n) - x(1:n - 1))*(y(2:n) + y(1:n - 1)))
    end function sum_trapez_1d

end module integrate

module integrate_substituted
    use integrate, only: integrand_i, integrate_1d
    use constants, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none

    private
    procedure(integrand_i), pointer, save :: integrand => null()
    real(dp), save :: left_limit, right_limit

    public :: integrate_1d_substituted

contains

    subroutine integrate_1d_substituted(f, a, b, result)
        procedure(integrand_i) :: f
        real(dp), intent(in) :: a, b
        real(dp), intent(out) :: result

        real(dp) :: sub_a, sub_b
        real(dp) :: left_result, right_result

        if (b < a) then
            print *, "Error in integrate_1d_substituted: b < a"
            error stop
        end if

        ! integration by globally adaptive interval subdivision (quadpack import)
        ! qags struggels with functions of form sqrt(x-a) and sqrt(b-x)
        ! Terefore we split the integral in the middle and substitute
        ! - left integral: t**2 = x - a
        ! - right integral: t**2 = b - x
        ! resulting in new limits
        ! - left integral: 0 to sqrt((b - a)*0.5_dp)
        ! - right integral: 0 to sqrt((b - a)*0.5_dp) (after absorbing [-] from trafo)
        ! the substitution handled with module states
        left_limit = a
        right_limit = b
        integrand => f
        sub_a = 0.0_dp
        sub_b = sqrt((b - a)*0.5_dp)
        call integrate_1d(left_substituted_integrand, sub_a, sub_b, left_result)
        call integrate_1d(right_substituted_integrand, sub_a, sub_b, right_result)
        left_limit = ieee_value(1.0_dp, ieee_quiet_nan)
        right_limit = ieee_value(1.0_dp, ieee_quiet_nan)
        integrand => null()

        result = left_result + right_result

    end subroutine integrate_1d_substituted

    function left_substituted_integrand(t)
        real(dp), intent(in) :: t
        real(dp) :: left_substituted_integrand

        real(dp) :: x

        x = t**2.0_dp + left_limit
        left_substituted_integrand = integrand(x)*2.0_dp*t
    end function left_substituted_integrand

    function right_substituted_integrand(t)
        real(dp), intent(in) :: t
        real(dp) :: right_substituted_integrand

        real(dp) :: x

        x = right_limit - t**2.0_dp
        right_substituted_integrand = integrand(x)*2.0_dp*t
    end function right_substituted_integrand

end module integrate_substituted
