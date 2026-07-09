program test_plot_anti_sigma
    use constants, only: dp, pi
    use utils, only: linspace
    use anti_sigma_field, only: anti_sigma_field_t
    use mock_perturbed_field, only: mock_perturbed_field_t
    use fieldline_mod, only: flock_of_fieldlines_t
    use make_fieldline, only: make_flock_from_labels

    use plot_quantities, only: plot_fieldlines_over_field

    implicit none

    real(dp), parameter :: M_pol = 0.0_dp, N_tor = 1.0_dp, nfp = N_tor
    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = -0.0125_dp, eps_1 = -0.0005_dp
    type(anti_sigma_field_t) :: field
    real(dp), parameter :: B_pert = -0.01_dp, M_pol_pert = 1.0_dp, N_tor_pert = 0.0_dp
    type(mock_perturbed_field_t) :: perturbed_field

    integer, parameter :: n_fieldlines = 20

    real(dp), dimension(n_fieldlines) :: theta_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    real(dp), parameter :: iota = 1.0_dp
    type(flock_of_fieldlines_t) :: flock

    call field%anti_sigma_field_init(M_pol, N_tor, B_0, eps_0, eps_1)
    call perturbed_field%mock_perturbed_field_init(field, &
                                                   M_pol_pert, &
                                                   N_tor_pert, &
                                                   B_pert)
    call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines + 1, temp)
    theta_0 = temp(1:n_fieldlines)

    call make_flock_from_labels(flock, &
                                theta_0, &
                                iota, &
                                perturbed_field, &
                                M_pol, &
                                N_tor, &
                                nfp)

    call plot_fieldlines_over_field(flock%fieldlines, perturbed_field)

end program test_plot_anti_sigma
