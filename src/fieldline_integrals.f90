module fieldline_integrals
    use constants, only: dp
    use fieldline_mod, only: fieldline_t
    use field_base, only: field_t
    implicit none

    private
    class(field_t), pointer, save :: current_field => null()
    type(fieldline_t), pointer, save :: current_fieldline => null()

    public :: calc_fieldline_integrals

contains

    subroutine calc_fieldline_integrals(field, fieldline)
        use integrate_substituted, only: integrate_1d_substituted
        type(fieldline_t), intent(inout), target :: fieldline
        class(field_t), intent(in), target :: field

        current_field => field
        current_fieldline => fieldline

        call integrate_1d_substituted(wrapper_local_radial_drift, &
                                      fieldline%phi_max(1), &
                                      fieldline%phi_max(2), &
                                      fieldline%radial_drift)

        call integrate_1d_substituted(wrapper_lambda_over_B_squared, &
                                      fieldline%phi_max(1), &
                                      fieldline%phi_max(2), &
                                      fieldline%integral_lambda_b_over_B_squared)

        call integrate_1d_substituted(wrapper_one_over_B_squared, &
                                      fieldline%phi_max(1), &
                                      fieldline%phi_max(2), &
                                      fieldline%integral_one_over_B_squared)

        call integrate_1d_substituted(wrapper_nabla_s_over_B_squared, &
                                      fieldline%phi_max(1), &
                                      fieldline%phi_max(2), &
                                      fieldline%integral_nabla_s_over_B_squared)

        current_field => null()
        current_fieldline => null()

    end subroutine calc_fieldline_integrals

    function wrapper_lambda_over_B_squared(phi)
        use fieldline_integrands, only: lambda_over_B_squared

        real(dp), intent(in) :: phi
        real(dp) :: wrapper_lambda_over_B_squared

        real(dp) :: theta

        theta = current_fieldline%get_theta(phi)
        wrapper_lambda_over_B_squared = lambda_over_B_squared(current_field, &
                                                              theta, &
                                                              phi, &
                                                              current_fieldline%eta_b)
    end function wrapper_lambda_over_B_squared

    function wrapper_local_radial_drift(phi)
        use fieldline_integrands, only: local_radial_drift

        real(dp), intent(in) :: phi
        real(dp) :: wrapper_local_radial_drift

        real(dp) :: theta

        theta = current_fieldline%get_theta(phi)
        wrapper_local_radial_drift = local_radial_drift(current_field, &
                                                        theta, &
                                                        phi, &
                                                        current_fieldline%eta_b)
    end function wrapper_local_radial_drift

    function wrapper_one_over_B_squared(phi)
        use fieldline_integrands, only: B_squared

        real(dp), intent(in) :: phi
        real(dp) :: wrapper_one_over_B_squared

        real(dp) :: theta

        theta = current_fieldline%get_theta(phi)
        wrapper_one_over_B_squared = 1.0_dp/B_squared(current_field, theta, phi)
    end function wrapper_one_over_B_squared

    function wrapper_nabla_s_over_B_squared(phi)
        use fieldline_integrands, only: nabla_s_over_B_squared

        real(dp), intent(in) :: phi
        real(dp) :: wrapper_nabla_s_over_B_squared

        real(dp) :: theta

        theta = current_fieldline%get_theta(phi)
        wrapper_nabla_s_over_B_squared = nabla_s_over_B_squared(current_field, &
                                                                theta, &
                                                                phi)
    end function wrapper_nabla_s_over_B_squared

end module fieldline_integrals
