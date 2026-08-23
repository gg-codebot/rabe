module shaing_callen_wrappers
    use constants, only: dp
    use fieldline_mod, only: fieldline_t
    use field_base, only: field_t
    implicit none
    private

    class(field_t), allocatable :: this_field
    type(fieldline_t) :: this_fieldline
    type(fieldline_t) :: null_fieldline
    real(dp) :: this_eta
    real(dp), parameter :: null_eta = -1.0_dp

    public :: this_field, this_fieldline, null_fieldline, this_eta, null_eta
    public :: wrapper_lambda_over_B_squared

contains

    function wrapper_lambda_over_B_squared(phi)
        use fieldline_integrands, only: lambda_over_B_squared

        real(dp), intent(in) :: phi
        real(dp) :: wrapper_lambda_over_B_squared

        real(dp) :: theta

        theta = this_fieldline%get_theta(phi)
        wrapper_lambda_over_B_squared = lambda_over_B_squared(this_field, &
                                                              theta, &
                                                              phi, &
                                                              this_eta)
    end function wrapper_lambda_over_B_squared

end module shaing_callen_wrappers
