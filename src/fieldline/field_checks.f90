module field_checks
    use constants, only: dp, pi
    use field_base, only: field_t

    implicit none
    private

    public :: suspect_omnigenous_origin_not_minimum

contains

    function suspect_omnigenous_origin_not_minimum(field, N_tor, M_pol, retol)
        use find_extrema, only: find_local_minima
        use find_extrema, only: find_global_extrema
        use field_base, only: field_t
        use fieldline_mod, only: fieldline_t
        use field_along_fieldline, only: set_field_and_fieldline
        use field_along_fieldline, only: unset_field_and_fieldline
        use field_along_fieldline, only: B_mod_along_fieldline, dB_dphi_along_fieldline
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: N_tor, M_pol
        real(dp), intent(in), optional :: retol
        logical :: suspect_omnigenous_origin_not_minimum

        real(dp) :: iota
        real(dp), dimension(2) :: chi_interval, phi_interval
        type(fieldline_t) :: fieldline
        real(dp), dimension(:), allocatable :: phi, phi_error, B, dB_dphi
        integer :: n_min, n_extrema, idx
        real(dp) :: phi_max, phi_tol
        real(dp), parameter :: height_retol = 0.3_dp
        real(dp), dimension(2) :: extrema
        real(dp) :: B_min, B_max, B_range, B_at_origin, height
        real(dp) :: B_min_error, B_max_error, B_range_error
        real(dp) :: height_error

        !> as M_pol and N_tor are whole numbers that must no both be zero
        !> M_pol*pi - N_tor should never zero
        if (abs(M_pol*pi - N_tor) < 1e-8) then
            print *, "Error: (M_pol*pi - N_tor) must not be (close) zero."
            print *, "abs(M_pol*iota - N_tor) = ", abs(M_pol*pi - N_tor)
            error stop
        end if

        iota = pi
        fieldline%iota = iota
        fieldline%phi_0 = 0.0_dp
        fieldline%theta_0 = 0.0_dp

        suspect_omnigenous_origin_not_minimum = .false.

        chi_interval = [-pi, pi]
        phi_interval = chi_interval/(M_pol*iota - N_tor)

        call set_field_and_fieldline(field, fieldline)

        call find_local_minima(B_mod_along_fieldline, phi_interval, phi, phi_error)
        n_min = size(phi)
        if (n_min == 0) then
            print *, "Warning in suspect_omnigenous_origin_not_minimum: ", &
                "No local dB_dphifound in well around origin! "
            suspect_omnigenous_origin_not_minimum = .true.
            call unset_field_and_fieldline()
            return
        end if
        allocate (B(n_min))
        allocate (dB_dphi(n_min))
        call B_mod_along_fieldline(phi, B)
        call dB_dphi_along_fieldline(phi, dB_dphi)
        idx = minloc(B, dim=1)
        B_min_error = abs(dB_dphi(idx))*phi_error(idx)
        B_min = B(idx)
        deallocate (B, dB_dphi)

        phi_tol = abs(phi_interval(2) - phi_interval(1))*1e-3_dp
        extrema = find_global_extrema(B_mod_along_fieldline, phi_interval, &
                                      abstol=phi_tol)
        n_extrema = size(extrema)
        allocate (B(n_extrema))
        allocate (dB_dphi(n_extrema))
        call B_mod_along_fieldline(extrema, B)
        call dB_dphi_along_fieldline(extrema, dB_dphi)
        idx = maxloc(B, dim=1)
        B_max_error = abs(dB_dphi(idx))*phi_tol
        B_max = B(idx)
        deallocate (B, dB_dphi)

        call unset_field_and_fieldline()

        !> is the B-difference of origin and minimum significant
        !> compared to the B-range?
        B_range = B_max - B_min
        B_range_error = B_min_error + B_max_error
        call field%compute_B_mod(0.0_dp, 0.0_dp, B_at_origin)
        height = B_at_origin - B_min
        height_error = B_min_error
        if (present(retol)) height_error = height_error + retol*B_at_origin
        if ((height - height_error) > height_retol*(B_range + B_range_error)) then
            print *, "Detected that B at origin of provided field is"
            print *, "significantly above the minimum B i.e. difference > "
            print *, height_retol*100.0_dp, "% of the total B range!"
            print *, "(B_at_origin - B_min) = ", height
            print *, "with estimated error = ", height_error
            print *, "B_max = ", B_max, " B_min = ", B_min
            print *, "(B_max-B_min) = ", B_range
            print *, "with estimated error = ", B_range_error
            suspect_omnigenous_origin_not_minimum = .true.
        end if

    end function suspect_omnigenous_origin_not_minimum

end module field_checks
