module deviation
    use constants, only: dp, pi
    use field_base, only: field_t
    use fieldline_mod, only: flock_of_fieldlines_t

    implicit none
    private

    public :: calc_deviation

contains

    subroutine calc_deviation(flock, deviation_A, deviation_B)
        use fieldline_labels, only: fieldline_modes_t
        use fieldline_labels, only: allocate_modes
        use fieldline_labels, only: fourier_transform_over_label
        use surface_average_mod, only: surface_average_t
        use surface_average_mod, only: calc_surface_averages
        use fit_functions, only: S_A, S_B
        use error_handling, only: failed_sanity_check

        type(flock_of_fieldlines_t), intent(in) :: flock
        real(dp), intent(out) :: deviation_A, deviation_B

        real(dp), parameter :: tol = 1e-12

        type(fieldline_modes_t) :: modes
        real(dp) :: iota_p, eta_b
        type(surface_average_t) :: average
        real(dp) :: symmetric_remainder
        real(dp) :: B_squared_sqrtg

        logical :: any_has_sin_part

        call fourier_transform_over_label(flock, modes)

        any_has_sin_part = .false.
        if (has_sin_modes(modes%delta_aspect_ratio)) then
            print *, "error: non-vanishing sin part of delta aspect ratio: "
            print *, "sin part: ", sum(abs(modes%delta_aspect_ratio%sin_coeffs))
            print *, "cos part: ", sum(abs(modes%delta_aspect_ratio%cos_coeffs))
            any_has_sin_part = .true.
        end if
        if (has_sin_modes(modes%delta_eta)) then
            print *, "error: non-vanishing sin part of delta eta: "
            print *, "sin part: ", sum(abs(modes%delta_eta%sin_coeffs))
            print *, "cos part: ", sum(abs(modes%delta_eta%cos_coeffs))
            any_has_sin_part = .true.
        end if
        if (any_has_sin_part) call failed_sanity_check()

        call calc_surface_averages(flock, average)

        iota_p = flock%iota_p
        eta_b = flock%eta_b

        deviation_A = pi*sum(modes%radial_drift%sin_coeffs* &
                             modes%delta_aspect_ratio%cos_coeffs* &
                             S_A(iota_p*modes%delta_aspect_ratio%mode_numbers))

        symmetric_remainder = pi*sum(modes%radial_drift%cos_coeffs* &
                                     modes%delta_aspect_ratio%sin_coeffs* &
                                     S_A(iota_p*modes%delta_aspect_ratio%mode_numbers))

        if (abs(symmetric_remainder/deviation_A) > tol) then
            print *, "warning: non-vanishing symmetric part of deviation A: "
            print *, "symmetric: ", symmetric_remainder
            print *, "antisymmetric: ", deviation_A
            print *, "ratio: ", symmetric_remainder/deviation_A
        end if

        deviation_A = deviation_A*average%B_squared/average%lambda_b* &
                      sqrt(eta_b)*sqrt(flock%I_ref)/average%normalization

        deviation_B = pi*sum(modes%radial_drift%sin_coeffs* &
                             modes%delta_eta%cos_coeffs* &
                             S_B(iota_p*modes%delta_eta%mode_numbers))

        symmetric_remainder = pi*sum(modes%radial_drift%cos_coeffs* &
                                     modes%delta_eta%sin_coeffs* &
                                     S_B(iota_p*modes%delta_eta%mode_numbers))

        if (abs(symmetric_remainder/deviation_B) > tol) then
            print *, "warning: non-vanishing symmetric part of deviation B: "
            print *, "symmetric: ", symmetric_remainder
            print *, "antisymmetric: ", deviation_B
            print *, "ratio: ", symmetric_remainder/deviation_B
        end if

        deviation_B = deviation_B*average%B_squared/average%lambda_b*0.5_dp/ &
                      average%normalization

    end subroutine calc_deviation

    function has_sin_modes(modes)
        use fieldline_labels, only: modes_t
        type(modes_t), intent(in) :: modes
        logical :: has_sin_modes

        real(dp), parameter :: tol = 1e-3, numerical_zero = 1e-8
        real(dp) :: sum_sin, sum_cos

        sum_sin = sum(abs(modes%sin_coeffs))
        sum_cos = sum(abs(modes%cos_coeffs))

        if (sum_cos > numerical_zero) then
            has_sin_modes = sum_sin/sum_cos > tol
        else
            has_sin_modes = sum_sin > numerical_zero
        end if
    end function has_sin_modes

end module deviation
