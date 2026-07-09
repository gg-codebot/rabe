program plot_deviation_poloidal_anti_sigma
    use constants, only: dp, pi
    use utils, only: linspace
    use anti_sigma_field, only: anti_sigma_field_t
    use mock_perturbed_field, only: mock_perturbed_field_t
    use fieldline_mod, only: flock_of_fieldlines_t
    use make_fieldline, only: make_flock_from_labels
    use deviation, only: calc_deviation
    use fit_functions, only: S_A, S_B

    use plot_quantities, only: plot_fieldlines_over_field
    use plot_quantities, only: plot_delta_eta
    use plot_quantities, only: plot_delta_A
    use plot_quantities, only: plot_deviation

    implicit none

    real(dp), parameter :: M_pol = 0.0_dp, N_tor = 1.0_dp, nfp = N_tor
    real(dp), parameter :: psi_edge = abs(-0.00785398_dp)/(2.0_dp*pi) !Tm^2
    real(dp), parameter :: R = 1.00_dp
    real(dp), parameter :: J_pol_over_N_tor = -5.0_dp*1e6
    real(dp), parameter :: I_tor = 0.0_dp

    real(dp), parameter :: ds_dr = 0.387524_dp*100.0_dp !1/m
    real(dp), parameter :: dr_dpsi = 1.0_dp/(ds_dr*psi_edge)

    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = -0.125, eps_1 = -0.05_dp
    real(dp), parameter :: delta_B_1 = 2.0_dp*1e-4
    real(dp), parameter :: M_pol_pert = 1.0_dp, N_tor_pert = 0.0_dp
    real(dp), parameter :: eps_ratio = eps_1/abs(eps_0)
    real(dp), parameter :: delta_A_1 = 0.25_dp*eps_ratio* &
                           (1.0_dp + 16.0_dp/3.0_dp*abs(eps_0))
    real(dp), parameter :: B_max = B_0*(1.0_dp + abs(eps_0)) + abs(delta_B_1)
    real(dp), parameter :: delta_eta_1 = -delta_B_1/B_max**2.0_dp
    type(anti_sigma_field_t) :: field_unperturbed
    type(mock_perturbed_field_t) :: field

    integer, parameter :: n_fieldlines = 41

    real(dp), dimension(n_fieldlines) :: xi_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    real(dp), parameter :: iota = 0.1618_dp
    type(flock_of_fieldlines_t) :: flock

    real(dp) :: deviation_A, deviation_B
    integer, parameter :: n_points = 100
    real(dp), dimension(n_points) :: nu_star
    real(dp) :: covariant_factor
    real(dp) :: off_factor_A, off_factor_B

    real(dp), parameter :: iota_p = sign(pi, N_tor)* &
                           (N_tor*iota + M_pol)/(N_tor - iota*M_pol)* &
                           nfp/(N_tor**2.0_dp + M_pol**2.0_dp)
    real(dp) :: off_factor_A_analytic, off_factor_B_analytic

    logical, parameter :: should_plot_others = .true.

    call field_unperturbed%anti_sigma_field_init(M_pol, N_tor, B_0, eps_0, eps_1)
    call field%mock_perturbed_field_init(field_unperturbed, &
                                         M_pol_pert, &
                                         N_tor_pert, &
                                         delta_B_1)
    call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines + 1, temp)
    xi_0 = temp(1:n_fieldlines)

    call make_flock_from_labels(flock, &
                                xi_0, &
                                iota, &
                                field, &
                                M_pol, &
                                N_tor, &
                                nfp)

    if (should_plot_others) then
        call plot_fieldlines_over_field(flock%fieldlines, field)
        call plot_delta_eta(flock%fieldlines, flock%iota_p, delta_eta_1)
        call plot_delta_A(flock%fieldlines, delta_A_1)
    end if

    call calc_deviation(flock, deviation_A, deviation_B)

    covariant_factor = -2.0_dp*1e-7*(J_pol_over_N_tor*abs(N_tor) + I_tor*iota)
    off_factor_A = deviation_A*dr_dpsi*sqrt(covariant_factor)*sqrt(0.5_dp*R*pi)
    off_factor_B = deviation_B*0.5*R*pi*dr_dpsi

    off_factor_A_analytic = -dr_dpsi*sqrt(2.0_dp*pi)*R*B_0*abs(eps_1)**2.0_dp &
                            *N_tor*S_A(iota_p)/(2.0_dp*abs(eps_0))**0.75_dp &
                            /(N_tor - iota*M_pol)**1.5_dp*0.5_dp &
                            *(1.0_dp + 6.0_dp*abs(eps_0))

    off_factor_B_analytic = eps_1*dr_dpsi*R*pi*delta_B_1*S_B(iota_p)*0.25_dp

    call plot_deviation(off_factor_A, &
                        off_factor_B, &
                        off_factor_A_analytic, &
                        off_factor_B_analytic)

end program plot_deviation_poloidal_anti_sigma
