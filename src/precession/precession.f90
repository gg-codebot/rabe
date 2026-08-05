module precession
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use constants, only: dp, pi, machine_eps
    use field_base, only: field_t
    use field_base, only: field_3D_t
    use fieldline_mod, only: fieldline_t, flock_of_fieldlines_t
    use grid_mod, only: integration_grid_t
    use grid_mod, only: fieldline_for_precession_t
    use grid_mod, only: set_integration_grids
    use grid_mod, only: compute_bounce_integrals
    use grid_mod, only: set_splines
    use interpolate, only: SplineData1D, construct_splines_1d, evaluate_splines_1d_many

    implicit none

    logical, parameter :: ignore_magnetic_drift = .true.

    type :: bounce_integrals_t
        real(dp) :: M_pol, N_tor, nfp, iota
        real(dp), dimension(:), allocatable :: t, eta
        real(dp), dimension(:), allocatable :: avg_I_j
        real(dp), dimension(:), allocatable :: avg_bounce_time_weighted
        real(dp), dimension(:), allocatable :: avg_magnetic_drift_weighted
        real(dp), dimension(:, :), allocatable :: radial_drift_sin
        integer :: n_modes
    end type bounce_integrals_t

contains

    subroutine setup_bounce_integrals(field, flock, s_tor, bounce_integrals, n_eta)
        use field_instance, only: initialize_field_instance
        use fourier, only: real_ft
        class(field_3D_t), intent(in) :: field
        type(flock_of_fieldlines_t), intent(in) :: flock
        real(dp), intent(in) :: s_tor
        type(bounce_integrals_t), intent(out) :: bounce_integrals
        integer, intent(in), optional :: n_eta

        type(fieldline_for_precession_t), dimension(:), allocatable :: fieldlines
        real(dp) :: eta_t, eta_c
        type(integration_grid_t) :: grid

        real(dp), dimension(:, :), allocatable :: radial_drift_weighted
        real(dp), dimension(:, :), allocatable :: radial_drift_cos

        integer :: n_fieldlines
        integer :: idx

        n_fieldlines = size(flock%fieldlines)
        allocate (fieldlines(n_fieldlines))
        fieldlines%fieldline_t = flock%fieldlines

        bounce_integrals%M_pol = flock%M_pol
        bounce_integrals%N_tor = flock%N_tor
        bounce_integrals%nfp = flock%nfp
        bounce_integrals%iota = flock%iota
        fieldlines%M_pol = bounce_integrals%M_pol
        fieldlines%N_tor = bounce_integrals%N_tor
        fieldlines%nfp = bounce_integrals%nfp
        if (ieee_is_nan(bounce_integrals%M_pol)) then
            print *, "M_pol: ", bounce_integrals%M_pol
            error stop "Error: M_pol is NaN."
        end if
        if (ieee_is_nan(bounce_integrals%N_tor)) then
            print *, "N_tor: ", bounce_integrals%N_tor
            error stop "Error: N_tor is NaN."
        end if
        if (ieee_is_nan(bounce_integrals%nfp)) then
            print *, "nfp: ", bounce_integrals%nfp
            error stop "Error: nfp is NaN."
        end if
        if (ieee_is_nan(bounce_integrals%iota)) then
            print *, "iota: ", bounce_integrals%iota
            error stop "Error: iota is NaN."
        end if
        if (abs(bounce_integrals%M_pol*bounce_integrals%iota - bounce_integrals%N_tor) &
            < machine_eps) then
            print *, "Error-Resonant: M_pol*iota - N_tor must not be zero."
            error stop
        end if

        eta_c = 1.0_dp/get_smallest_maximum(fieldlines)
        call set_fieldline_minima(field, fieldlines)
        eta_t = 1.0_dp/get_biggest_minimum(fieldlines)
        call set_integration_grids(eta_t, eta_c, grid, n_eta=n_eta)
        fieldlines%grid = grid

        call initialize_field_instance(field)
        allocate (radial_drift_weighted(n_fieldlines, grid%n))
        allocate (bounce_integrals%avg_bounce_time_weighted(grid%n))
        allocate (bounce_integrals%avg_I_j(grid%n))
        allocate (bounce_integrals%avg_magnetic_drift_weighted(grid%n))
        bounce_integrals%avg_bounce_time_weighted = 0.0_dp
        bounce_integrals%avg_I_j = 0.0_dp
        bounce_integrals%avg_magnetic_drift_weighted = 0.0_dp
        do idx = 1, n_fieldlines
            call compute_bounce_integrals(field, &
                                          fieldlines(idx), &
                                          s_tor, &
                                          fieldlines(idx)%grid)
            radial_drift_weighted(idx, :) = fieldlines(idx)%grid%radial_drift_weighted
            bounce_integrals%avg_bounce_time_weighted = &
                bounce_integrals%avg_bounce_time_weighted + &
                fieldlines(idx)%grid%bounce_time_weighted
            bounce_integrals%avg_I_j = bounce_integrals%avg_I_j + &
                                       fieldlines(idx)%grid%I_j
            bounce_integrals%avg_magnetic_drift_weighted = &
                bounce_integrals%avg_magnetic_drift_weighted + &
                fieldlines(idx)%grid%poloidal_drift_weighted
        end do
        if (ignore_magnetic_drift) then
            bounce_integrals%avg_magnetic_drift_weighted = 0.0_dp
        else
            bounce_integrals%avg_magnetic_drift_weighted = &
                bounce_integrals%avg_magnetic_drift_weighted/n_fieldlines
        end if
        bounce_integrals%avg_bounce_time_weighted = &
            bounce_integrals%avg_bounce_time_weighted/n_fieldlines
        bounce_integrals%avg_I_j = bounce_integrals%avg_I_j/n_fieldlines

        bounce_integrals%n_modes = n_fieldlines/2 + 1
        allocate (radial_drift_cos(bounce_integrals%n_modes, grid%n))
        allocate (bounce_integrals%radial_drift_sin(bounce_integrals%n_modes, grid%n))
        do idx = 1, grid%n
            call real_ft(fieldlines%xi_0, &
                         radial_drift_weighted(:, idx), &
                         radial_drift_cos(:, idx), &
                         bounce_integrals%radial_drift_sin(:, idx))
        end do

        bounce_integrals%t = grid%t
        bounce_integrals%eta = grid%eta

        deallocate (fieldlines)
        deallocate (radial_drift_weighted)
        deallocate (radial_drift_cos)

    end subroutine setup_bounce_integrals

    subroutine compute_flux_modes(bounce_integrals, l_c, Omega_hat, flux_modes)
        use splines_instance, only: initialize_splines
        use splines_instance, only: initialize_prefactor
        use splines_instance, only: initialize_radial_drift_spline
        use splines_instance, only: get_flux_mode
        use splines_instance, only: initialize_startup
        type(bounce_integrals_t), intent(in) :: bounce_integrals
        real(dp), intent(in) :: l_c
        real(dp), intent(in) :: Omega_hat
        real(dp), dimension(:), allocatable, intent(out) :: flux_modes

        real(dp), dimension(:), allocatable :: poloidal_drift_weighted
        real(dp), dimension(:), allocatable :: electric_drift_weighted
        real(dp) :: mode_factor
        integer :: idx

        allocate (electric_drift_weighted(size(bounce_integrals%t)))
        electric_drift_weighted = Omega_hat*bounce_integrals%avg_bounce_time_weighted
        poloidal_drift_weighted = bounce_integrals%avg_magnetic_drift_weighted + &
                                  electric_drift_weighted

        call initialize_splines(bounce_integrals%t, &
                                bounce_integrals%eta, &
                                bounce_integrals%avg_I_j, &
                                poloidal_drift_weighted)

        allocate (flux_modes(bounce_integrals%n_modes))
        flux_modes(1) = 0.0_dp
        do idx = 2, bounce_integrals%n_modes
            mode_factor = 0.5_dp*real(idx - 1, dp)*l_c*bounce_integrals%nfp/ &
                 (bounce_integrals%M_pol*bounce_integrals%iota - bounce_integrals%N_tor)
            call initialize_prefactor(mode_factor)
            call initialize_startup(bounce_integrals%t, bounce_integrals%eta, &
                                    bounce_integrals%avg_I_j, mode_factor)
            call initialize_radial_drift_spline(bounce_integrals%t, &
                                              bounce_integrals%radial_drift_sin(idx, :))
            call get_flux_mode(bounce_integrals%t(1), &
                          bounce_integrals%t(size(bounce_integrals%t)), flux_modes(idx))
        end do

        deallocate (poloidal_drift_weighted)
        deallocate (electric_drift_weighted)

    end subroutine compute_flux_modes

    function get_smallest_maximum(fieldlines) result(smallest_maximum)
        class(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp) :: smallest_maximum

        real(dp), dimension(size(fieldlines)) :: smallest_maximum_per_line
        integer :: idx

        do idx = 1, size(fieldlines)
            smallest_maximum_per_line(idx) = minval(fieldlines(idx)%B_max)
        end do
        smallest_maximum = minval(smallest_maximum_per_line)

        if (smallest_maximum <= 0.0_dp) then
            error stop "Error: Smallest maximum is not positiv."
        end if
        if (ieee_is_nan(smallest_maximum)) then
            error stop "Error: Smallest maximum is NaN."
        end if

    end function get_smallest_maximum

    subroutine set_fieldline_minima(field, fieldlines)
        class(field_t), intent(in) :: field
        class(fieldline_for_precession_t), dimension(:), intent(inout) :: fieldlines

        integer :: idx

        do idx = 1, size(fieldlines)
            call find_magnetic_well_bottom(field, &
                                           fieldlines(idx)%fieldline_t, &
                                           fieldlines(idx)%phi_min, &
                                           fieldlines(idx)%B_min)
        end do
    end subroutine set_fieldline_minima

    function get_biggest_minimum(fieldlines) result(biggest_minimum)
        class(fieldline_for_precession_t), dimension(:), intent(in) :: fieldlines
        real(dp) :: biggest_minimum

        biggest_minimum = maxval(fieldlines%B_min)

        if (biggest_minimum <= 0.0_dp) then
            error stop "Error: Biggest minimum is not positiv."
        end if
        if (ieee_is_nan(biggest_minimum)) then
            error stop "Error: Biggest minimum is NaN."
        end if
    end function get_biggest_minimum

    subroutine integrate_radial_drift(grid, fieldline, integral)
        use odeint_allroutines_sub, only: odeint_allroutines
        use grid_instance, only: initialize_grid_instance
        type(integration_grid_t), intent(in) :: grid
        type(fieldline_t), intent(in) :: fieldline
        real(dp), intent(out) :: integral

        integer, parameter :: ndim = 1
        real(dp) :: t_start, t_end
        real(dp), parameter :: relerr = 1e-8_dp
        real(dp), dimension(ndim) :: y

        t_start = grid%t(1)
        t_end = grid%t(size(grid%t))

        call initialize_grid_instance(grid)
        y(1) = 0.0_dp
        call odeint_allroutines(y, ndim, t_start, t_end, relerr, rhs)

        integral = y(1)
    end subroutine integrate_radial_drift

    subroutine rhs(t, y, dydt)
        use grid_instance, only: get_radial_drift
        real(dp), intent(in) :: t
        real(dp), dimension(:), intent(in) :: y
        real(dp), dimension(:), intent(out) :: dydt

        real(dp) :: radial_drift
        call get_radial_drift(t, dydt(1))

    end subroutine rhs

    function get_fieldline_at_global_maximum(fieldlines) result(max_fieldline)
        class(fieldline_t), dimension(:), intent(in) :: fieldlines
        type(fieldline_t) :: max_fieldline

        integer :: max_index

        max_index = maxloc(fieldlines%B_max(1), dim=1)
        max_fieldline = fieldlines(max_index)
    end function get_fieldline_at_global_maximum

    subroutine find_magnetic_well_bottom(field, fieldline, phi_bottom, B_bottom)
        use find_extrema, only: find_local_minima
        class(field_t), intent(in) :: field
        class(fieldline_t), intent(in) :: fieldline
        real(dp), intent(out) :: phi_bottom, B_bottom

        integer :: found_minima, of_smallest_B
        real(dp), dimension(:), allocatable :: phi_min, B_min

        call find_local_minima(B_along_fieldline, fieldline%phi_max, phi_min)

        found_minima = size(phi_min)
        allocate (B_min(found_minima))
        if (found_minima < 1) then
            print *, "---------------------------------------------------------"
            print *, "---------------------------------------------------------"
            print *, "---------------------------------------------------------"
            print *, "error in find_magnetic_well_bottom: "
            print *, "Found no local minimum in provided interval!"
            print *, "theta_0: ", fieldline%theta_0
            print *, "phi_0: ", fieldline%phi_0
            print *, "interval: ", fieldline%phi_max
            print *, "---------------------------------------------------------"
            print *, "---------------------------------------------------------"
            print *, "---------------------------------------------------------"
            error stop
        elseif (found_minima > 2) then
            print *, "---------------------------------------------------------"
            print *, "---------------------------------------------------------"
            print *, "---------------------------------------------------------"
            print *, "warning in find_magnetic_well_bottom: "
            print *, "More than one local minimum in provided interval!"
            print *, "theta_0: ", fieldline%theta_0
            print *, "phi_0: ", fieldline%phi_0
            print *, "interval: ", fieldline%phi_max
            print *, "phi_min: ", phi_min(1:found_minima)
            print *, "---------------------------------------------------------"
            print *, "---------------------------------------------------------"
            print *, "---------------------------------------------------------"
        end if

        call B_along_fieldline(phi_min(1:found_minima), B_min(1:found_minima))
        of_smallest_B = minloc(B_min(1:found_minima), dim=1)
        phi_bottom = phi_min(of_smallest_B)
        B_bottom = B_min(of_smallest_B)

    contains

        subroutine B_along_fieldline(phi, B)
            real(dp), dimension(:), intent(in) :: phi
            real(dp), dimension(:), intent(out) :: B

            real(dp), dimension(size(phi)) :: theta
            integer :: idx

            theta = fieldline%get_theta(phi)
            do idx = 1, size(phi)
                call field%compute_B_mod(theta(idx), phi(idx), B(idx))
            end do
        end subroutine B_along_fieldline

    end subroutine find_magnetic_well_bottom

    function find_turning_points(field, fieldline, eta, reltol_in) result(phi_turning)
        use find_extrema, only: find_global_minimum
        class(field_t), intent(in) :: field
        type(fieldline_for_precession_t), intent(in) :: fieldline
        real(dp), intent(in) :: eta
        real(dp), intent(in), optional :: reltol_in
        real(dp), dimension(2) :: phi_turning, extremum_locations
        real(dp) :: interval(2)

        real(dp) :: reltol

        if (present(reltol_in)) then
            reltol = reltol_in
        else
            reltol = 1e-4_dp
        end if

        ! turning points are shifted slightly away from the extrema to avoid
        ! numerical issues with the integrand
        interval(1) = fieldline%phi_max(1)
        interval(2) = fieldline%phi_min
        extremum_locations = find_global_minimum(lambda_squared_along_fieldline, interval, reltol)
        phi_turning(1) = extremum_locations(1) + (interval(2) - interval(1))*reltol

        interval(1) = fieldline%phi_min
        interval(2) = fieldline%phi_max(2)
        extremum_locations = find_global_minimum(lambda_squared_along_fieldline, interval, reltol)
        phi_turning(2) = extremum_locations(1) - (interval(2) - interval(1))*reltol

    contains

        subroutine lambda_squared_along_fieldline(phi, lambda_squared)
            real(dp), dimension(:), intent(in) :: phi
            real(dp), dimension(:), intent(out) :: lambda_squared

            real(dp), dimension(size(phi)) :: theta, B
            integer :: idx

            theta = fieldline%get_theta(phi)
            do idx = 1, size(phi)
                call field%compute_B_mod(theta(idx), phi(idx), B(idx))
            end do
            lambda_squared = abs(1.0_dp - B*eta)
        end subroutine lambda_squared_along_fieldline

    end function find_turning_points

end module precession
