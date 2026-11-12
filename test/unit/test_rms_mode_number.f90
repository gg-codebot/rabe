program test_rms_mode_number
    use constants, only: dp
    use fieldline_labels, only: modes_t, allocate_modes, calc_rms_mode_number

    implicit none

    real(dp), parameter :: reltol = 1e-14
    type(modes_t) :: modes
    real(dp) :: expected

    ! a single mode gives its mode number, for cos and sin alike
    call allocate_modes(modes, 6)
    modes%sin_coeffs(3) = -0.7_dp
    call check(calc_rms_mode_number(modes), 2.0_dp, "single sin mode m=2")
    modes%sin_coeffs(3) = 0.0_dp
    modes%cos_coeffs(5) = 1.3_dp
    call check(calc_rms_mode_number(modes), 4.0_dp, "single cos mode m=4")

    ! two modes: power weighted, independent of the phases
    modes%cos_coeffs = 0.0_dp
    modes%sin_coeffs(2) = 1.0_dp
    modes%cos_coeffs(4) = -0.5_dp
    expected = sqrt((1.0_dp*1.0_dp + 9.0_dp*0.25_dp)/(1.0_dp + 0.25_dp))
    call check(calc_rms_mode_number(modes), expected, "modes m=1 and m=3")

contains

    subroutine check(found, wanted, case_name)
        real(dp), intent(in) :: found, wanted
        character(len=*), intent(in) :: case_name

        if (abs(found - wanted) > reltol*abs(wanted)) then
            print *, "test_rms_mode_number failed for ", case_name
            print *, "found: ", found, " expected: ", wanted
            error stop
        end if
    end subroutine check
end program test_rms_mode_number
