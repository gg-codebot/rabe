program plot_integrate_radial_drift
    use constants, only: dp, pi
    use boozer_field, only: boozer_field_t
    use fieldline_mod, only: flock_of_fieldlines_t
    use make_fieldline, only: make_flock_from_labels
    use grid_mod, only: integration_grid_t
    use grid_mod, only: fieldline_for_precession_t
    use grid_mod, only: set_integration_grids
    use grid_mod, only: compute_bounce_integrals
    use grid_mod, only: set_splines
    use precession, only: get_fieldline_at_global_maximum
    use precession, only: find_magnetic_well_bottom
    use precession, only: integrate_radial_drift
    use field_instance, only: initialize_field_instance
    use utils, only: linspace

    implicit none

    real(dp), parameter :: M_pol = 0.0_dp, N_tor = 4.0_dp
    integer, parameter :: n_fieldlines = 50
    real(dp), parameter :: s_tor = 0.25_dp
    real(dp), parameter :: nfp = N_tor

    real(dp), dimension(n_fieldlines) :: xi_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    real(dp) :: iota
    real(dp) :: B_covariant_phi, B_covariant_theta
    real(dp) :: phi_bottom, B_min
    real(dp) :: lowest_B_max, eta_t, eta_c
    real(dp) :: radial_drift_integral
    real(dp) :: check

    type(flock_of_fieldlines_t) :: flock
    type(fieldline_for_precession_t) :: precession_fieldline
    type(integration_grid_t) :: grid
    type(boozer_field_t) :: field

    character(len=*), parameter :: nc_file = "input/wout_squid_20230921_v1_shifted.nc"

    call field%init_from_vmec(nc_file, grid_refinement=5)
    call field%fix_to_surface(s_tor)
    call field%get_iota(s_tor, iota)
    call field%get_covariant_components(B_covariant_theta, B_covariant_phi)

    call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines + 1, temp)
    xi_0 = temp(1:n_fieldlines)

    call make_flock_from_labels(flock, &
                                xi_0, &
                                iota, &
                                field, &
                                M_pol, &
                                N_tor, &
                                nfp)

    precession_fieldline%fieldline_t = get_fieldline_at_global_maximum(flock%fieldlines)
    precession_fieldline%M_pol = M_pol
    precession_fieldline%N_tor = N_tor
    precession_fieldline%nfp = nfp

    call find_magnetic_well_bottom(field, precession_fieldline, phi_bottom, B_min)
    lowest_B_max = minval(precession_fieldline%B_max)
    eta_t = 1.0_dp/B_min
    eta_c = 1.0_dp/lowest_B_max
    precession_fieldline%phi_min = phi_bottom
    precession_fieldline%B_min = B_min

    call set_integration_grids(eta_t, eta_c, grid)
    call initialize_field_instance(field)
    call compute_bounce_integrals(field, precession_fieldline, s_tor, grid)
    call set_splines(grid)

    call integrate_radial_drift(grid, precession_fieldline%fieldline_t, radial_drift_integral)
    print *, "integrate_radial_drift: ", radial_drift_integral
    check = precession_fieldline%radial_drift
    check = check/0.75_dp*(B_covariant_phi + iota*B_covariant_theta)/field%psi_tor_edge
    print *, "check: ", check

end program plot_integrate_radial_drift
