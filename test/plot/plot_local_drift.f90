program plot_anti_sigma
    use constants, only: dp, pi
    use utils, only: linspace
    use anti_sigma_field, only: anti_sigma_field_t
    use mock_perturbed_field, only: mock_perturbed_field_t
    use fieldline_mod, only: fieldline_t

    use plot_quantities, only: plot_fieldline_over_local_drift
    use plot_quantities, only: plot_local_drift_over_fieldline

    implicit none

    real(dp), parameter :: M_pol = -1.0_dp, N_tor = 1.0_dp, nfp = N_tor
    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = -0.0125_dp, eps_1 = -0.0005_dp
    type(anti_sigma_field_t) :: field
    real(dp), parameter :: B_pert = -0.01_dp, M_pol_pert = 1.0_dp, N_tor_pert = 0.0_dp
    type(mock_perturbed_field_t) :: perturbed_field
    real(dp), parameter :: B_max = B_0*(1.0_dp + abs(eps_0) + abs(eps_1)) + abs(B_pert)
    real(dp), parameter :: B_min = B_0*(1.0_dp - abs(eps_0) - abs(eps_1)) - abs(B_pert)

    real(dp), parameter :: iota = 1.0_dp
    integer, parameter :: n_eta = 100
    real(dp), dimension(n_eta) :: etas
    real(dp), dimension(2), parameter :: interval = [0.0_dp, 2.0_dp*pi/nfp]
    type(fieldline_t) :: fieldline

    call field%anti_sigma_field_init(M_pol, N_tor, B_0, eps_0, eps_1)
    call perturbed_field%mock_perturbed_field_init(field, &
                                                   M_pol_pert, &
                                                   N_tor_pert, &
                                                   B_pert)
    fieldline%iota = iota
    fieldline%phi_0 = pi/nfp
    fieldline%theta_0 = 0.0_dp
    fieldline%xi_0 = pi*nfp*M_pol/(M_pol**2 + N_tor**2)

    call linspace(1.0_dp/B_max, 1.0_dp/B_min, n_eta, etas)
    call plot_fieldline_over_local_drift(fieldline, perturbed_field, etas(1), interval)
    call plot_local_drift_over_fieldline(perturbed_field, fieldline, etas, &
                                         interval=interval)

end program plot_anti_sigma
