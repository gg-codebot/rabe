program plot_bounce_time_averages
    use constants, only: dp, pi, machine_eps
    use anti_sigma_field, only: anti_sigma_field_t
    use mock_perturbed_field, only: mock_perturbed_field_t
    use mock_field_3d, only: mock_field_3d_t
    use boozer_field, only: boozer_field_t
    use fieldline_mod, only: flock_of_fieldlines_t
    use make_fieldline, only: make_flock_from_labels
    use grid_mod, only: integration_grid_t
    use grid_mod, only: fieldline_for_precession_t
    use precession_reference, only: get_fieldline_at_global_maximum
    use precession, only: find_magnetic_well_bottom
    use grid_mod, only: set_integration_grids
    use grid_mod, only: compute_bounce_integrals
    use grid_mod, only: set_splines
    use precession_reference, only: evaluate_grid_splines
    use field_instance, only: initialize_field_instance
    use utils, only: linspace
    use myplot_module, only: myplot

    use plot_quantities, only: plot_local_drift_over_fieldline

    implicit none

    real(dp), parameter :: M_pol = 0.0_dp, N_tor = 4.0_dp
    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = -0.0125_dp, eps_1 = -0.0005_dp
    type(anti_sigma_field_t) :: field_2D
    real(dp), parameter :: B_pert = 0.001_dp, M_pol_pert = 1.0_dp, N_tor_pert = 0.0_dp
    type(mock_perturbed_field_t) :: perturbed_field_2D
    type(mock_field_3d_t) :: field_mock

    integer, parameter :: n_fieldlines = 50
    real(dp), dimension(n_fieldlines) :: xi_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    real(dp) :: iota
    real(dp), parameter :: nfp = N_tor
    type(flock_of_fieldlines_t) :: flock
    type(fieldline_for_precession_t) :: precession_fieldline
    type(integration_grid_t) :: grid

    real(dp), parameter :: s_tor = 0.25_dp
    integer, parameter :: n_spline_plot = 300
    real(dp), dimension(:), allocatable :: eta
    real(dp), dimension(:), allocatable :: bounce_time
    real(dp), dimension(:), allocatable :: bounce_time_deep
    real(dp), dimension(:), allocatable :: I_j, I_j_boundary
    real(dp), dimension(:), allocatable :: t_spline
    real(dp), dimension(:), allocatable :: eta_spline
    real(dp), dimension(:), allocatable :: I_j_spline_vals
    real(dp), dimension(:), allocatable :: bounce_coef_spline_vals
    real(dp), dimension(:), allocatable :: radial_drift_spline_vals
    real(dp), dimension(:), allocatable :: poloidal_drift_spline_vals
    real(dp) :: B_min, eta_b, phi_bottom
    real(dp) :: lowest_B_max, eta_t, eta_c, theta_min
    real(dp) :: well_depth
    real(dp) :: dummy(2), average_B, sqrtg
    real(dp), dimension(3) :: x, b_der, h_covar, h_ctrvr, h_curl
    integer, dimension(2), parameter :: figsize = [14, 10]

    character(len=*), parameter :: nc_file = "input/wout_squid_20230921_v1_shifted.nc"
    type(boozer_field_t) :: field

    type(myplot) :: plt

    ! call field_2D%anti_sigma_field_init(M_pol, N_tor, B_0, eps_0, eps_1)
    ! call perturbed_field_2D%mock_perturbed_field_init(field_2D, &
    !                                                   M_pol_pert, &
    !                                                   N_tor_pert, &
    !                                                   B_pert)
    ! call field%mock_field_3d_init(perturbed_field_2D)
    ! call field%get_iota(s_tor, iota)

    call field%init_from_vmec(nc_file, grid_refinement=5)
    call field%fix_to_surface(s_tor)
    call field%get_iota(s_tor, iota)

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
    eta_b = precession_fieldline%eta_b

    call set_integration_grids(eta_t, eta_c, grid)

    call plot_local_drift_over_fieldline(field, precession_fieldline%fieldline_t, grid%eta, precession_fieldline%phi_max)

    call initialize_field_instance(field)
    call compute_bounce_integrals(field, precession_fieldline, s_tor, grid)
    call set_splines(grid)

    allocate (eta(grid%n - 1))
    allocate (bounce_time(size(eta)))
    allocate (I_j(size(eta)))

    eta = grid%eta(2:grid%n)

    bounce_time = grid%normalized_bounce_time(2:grid%n)

    I_j = grid%I_j(2:grid%n)

    theta_min = precession_fieldline%get_theta(phi_bottom)
    x = [s_tor, theta_min, phi_bottom]
    call field%evaluate(x, dummy(1), sqrtg, b_der, h_covar, h_ctrvr, h_curl)

    allocate (bounce_time_deep(size(eta)))
    allocate (I_j_boundary(size(eta)))
    allocate (t_spline(n_spline_plot))
    allocate (eta_spline(n_spline_plot))
    allocate (I_j_spline_vals(n_spline_plot))
    allocate (bounce_coef_spline_vals(n_spline_plot))
    allocate (radial_drift_spline_vals(n_spline_plot))
    allocate (poloidal_drift_spline_vals(n_spline_plot))

    well_depth = max(machine_eps, 1.0_dp - B_min*eta_b)
    average_B = 0.5_dp*(B_min + lowest_B_max)
    bounce_time_deep = 4.0_dp*pi/abs(h_ctrvr(2))/sqrt(well_depth)/ &
                       abs(M_pol*iota - N_tor)
    I_j_boundary = 4.0_dp*sqrt(well_depth)/(average_B*h_ctrvr(3))/ &
                   abs(M_pol*iota - N_tor)

    call linspace(grid%t(1), grid%t(grid%n), n_spline_plot, t_spline)
    eta_spline = eta_c + t_spline**2.0_dp
    call evaluate_grid_splines(grid, t_spline, I_j_spline_vals, &
                               bounce_coef_spline_vals, &
                               radial_drift_spline_vals, poloidal_drift_spline_vals)

    call plt%initialize(xlabel="$1 - \eta/\eta_t$", &
                        ylabel="$\tau_b v_{\mathrm{th}}$", &
                        legend=.true., &
                        figsize=figsize, &
                        title="$\tau_b$ vs $\eta/\eta_c - 1$")

    call plt%add_plot(eta/eta_c - 1.0_dp, bounce_time, &
                      label="$\tau_b$ (numerical)", &
                      linestyle="r-")
    call plt%add_plot(eta/eta_c - 1.0_dp, bounce_time_deep, &
                      label="$\tau_b$ (deep-trapped estimate)", &
                      linestyle="k--")
    call plt%show()

    call plt%initialize(xlabel="$\eta/\eta_c - 1$", &
                        ylabel="$I_j$", &
                        legend=.true., &
                        figsize=figsize, &
                        title="$I_j$ vs $\eta/\eta_c - 1$")

    call plt%add_plot(eta/eta_c - 1.0_dp, I_j, &
                      label="$I_j$ (numerical)", &
                      linestyle="r-")
    call plt%add_plot(eta/eta_c - 1.0_dp, I_j_boundary, &
                      label="$I_j$ (trapped-passing boundary estimate)", &
                      linestyle="k--")
    call plt%add_plot(eta_spline/eta_c - 1.0_dp, I_j_spline_vals, &
                      label="$I_j$ (spline)", &
                      linestyle="b--")

    call plt%show()

    call plt%initialize(xlabel="$t$", &
                        ylabel="$C_{\mathrm{bounce}}$", &
                        legend=.true., &
                        figsize=figsize, &
                        title="$C_{\mathrm{bounce}}$ vs $t$")

    call plt%add_plot(grid%t, grid%bounce_time_weighted, &
                      label="$C_{\mathrm{bounce}}$ (numerical)", &
                      linestyle="r-")
    call plt%add_plot(t_spline, bounce_coef_spline_vals, &
                      label="$C_{\mathrm{bounce}}$ (spline)", &
                      linestyle="b--")

    call plt%show()

    call plt%initialize(xlabel="$t$", &
                        ylabel="$D_{r}$", &
                        legend=.true., &
                        figsize=figsize, &
                        title="$D_{r}$ vs $t$")

    call plt%add_plot(grid%t(2:grid%n), &
                      grid%normalized_radial_drift(2:grid%n), &
                      label="$D_{r}$ (numerical)", &
                      linestyle="r-")

    call plt%show()

    call plt%initialize(xlabel="$t$", &
                        ylabel="$D_{r}$", &
                        legend=.true., &
                        figsize=figsize, &
                        title="$D_{\vartheta}$ vs $t$")

    call plt%add_plot(grid%t(2:grid%n), grid%poloidal_drift(2:grid%n), &
                      label="$D_{\vartheta}$ (numerical)", &
                      linestyle="r-")

    call plt%show()

    call plt%initialize(xlabel="$t$", &
                        ylabel="$D_{r,\mathrm{w}}$", &
                        legend=.true., &
                        figsize=figsize, &
                        title="$D_{r,\mathrm{w}}$ vs $t$")

    call plt%add_plot(grid%t, grid%radial_drift_weighted, &
                      label="$D_{r,\mathrm{w}}$ (numerical)", &
                      linestyle="r-")
    call plt%add_plot(t_spline, radial_drift_spline_vals, &
                      label="$D_{r,\mathrm{w}}$ (spline)", &
                      linestyle="b--")

    call plt%show()

    call plt%initialize(xlabel="$t$", &
                        ylabel="$D_{\vartheta,\mathrm{w}}$", &
                        legend=.true., &
                        figsize=figsize, &
                        title="$D_{\vartheta,\mathrm{w}}$ vs $t$")

    call plt%add_plot(grid%t, grid%poloidal_drift_weighted, &
                      label="$D_{\vartheta,\mathrm{w}}$ (numerical)", &
                      linestyle="r-")
    call plt%add_plot(t_spline, poloidal_drift_spline_vals, &
                      label="$D_{\vartheta,\mathrm{w}}$ (spline)", &
                      linestyle="b--")

    call plt%show()

end program plot_bounce_time_averages
