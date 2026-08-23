module fieldline_integrands
    use constants, only: dp
    use field_base, only: field_t

    implicit none
    private

    public :: local_radial_drift, lambda_over_B_squared, B_squared
    public :: calc_lambda_squared, nabla_s_over_B_squared

contains

    function local_radial_drift(field, theta, phi, eta)
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: theta, phi, eta
        real(dp) :: local_radial_drift

        real(dp) :: B, dB_dx(3), lambda_squared, dB_dtheta

        call field%compute_B_and_dB_dx(theta, phi, B, dB_dx)
        lambda_squared = calc_lambda_squared(B, eta)

        dB_dtheta = dB_dx(2)

        local_radial_drift = 0.5_dp*sqrt(lambda_squared)/B**3.0_dp &
                             *(3.0_dp + lambda_squared) &
                             *dB_dtheta
    end function local_radial_drift

    function lambda_over_B_squared(field, theta, phi, eta)
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: theta, phi, eta
        real(dp) :: lambda_over_B_squared

        real(dp) :: B

        call field%compute_B_mod(theta, phi, B)
        lambda_over_B_squared = sqrt(calc_lambda_squared(B, eta))/B**2.0_dp
    end function lambda_over_B_squared

    function B_squared(field, theta, phi)
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: theta, phi
        real(dp) :: B_squared

        real(dp) :: B

        call field%compute_B_mod(theta, phi, B)
        B_squared = B**2.0_dp
    end function B_squared

    function calc_lambda_squared(B, eta) result(lambda_squared)
        real(dp), intent(in) :: B, eta
        real(dp) :: lambda_squared

        lambda_squared = 1.0_dp - B*eta

        if (lambda_squared .lt. 0.0_dp) then
            print *, "Square of pitch parameter (1 - B*eta) is negative!"
            print *, "B: ", B
            print *, "eta: ", eta
            print *, "(1 - B*eta): ", lambda_squared
            error stop
        end if
    end function calc_lambda_squared

    function nabla_s_over_B_squared(field, theta, phi)
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: theta, phi
        real(dp) :: nabla_s_over_B_squared

        real(dp) :: B, nabla_s

        call field%compute_B_mod(theta, phi, B)
        call field%compute_nabla_s(theta, phi, nabla_s)
        nabla_s_over_B_squared = nabla_s/B**2.0_dp
    end function nabla_s_over_B_squared

end module fieldline_integrands
