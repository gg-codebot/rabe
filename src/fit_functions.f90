module fit_functions
    use constants, only: dp, pi
    implicit none
    private

    real(dp), parameter :: eps = 1e-15_dp

    public :: S_A, S_B

contains

    elemental function S_A(angle)
        real(dp), intent(in) :: angle
        real(dp) :: S_A

        real(dp) :: angle_in_period
        real(dp), parameter :: fit_fac = 0.3_dp

        angle_in_period = modulo(angle + pi, 2.0_dp*pi) - pi

        S_A = fit_fac*(angle_in_period - sign(0.5_dp*pi, angle_in_period))
        if (abs(mod(angle, pi)) < eps) S_A = 0.0_dp
    end function S_A

    elemental function S_B(angle)
        real(dp), intent(in) :: angle
        real(dp) :: S_B

        real(dp) :: angle_in_period
        real(dp), parameter :: fit_fac = 2.0_dp

        angle_in_period = modulo(angle + pi, 2.0_dp*pi) - pi

        S_B = sign(fit_fac, angle_in_period)
        if (abs(mod(angle, pi)) < eps) S_B = 0.0_dp
    end function S_B

end module fit_functions
