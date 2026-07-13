program test_precession_analytic
    use constants, only: dp, pi
    use anti_sigma_field, only: anti_sigma_field_t
    use mock_perturbed_field, only: mock_perturbed_field_t
    use mock_field_3d, only: mock_field_3d_t
    use fieldline_mod, only: flock_of_fieldlines_t
    use make_fieldline, only: make_flock_from_labels
    use precession, only: bounce_integrals_t
    use precession, only: setup_bounce_integrals
    use precession, only: compute_flux_modes
    use utils, only: linspace

    implicit none

    real(dp), parameter :: M_pol = 0.0_dp, N_tor = 1.0_dp
    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = -0.0125_dp, &
                           eps_1 = -0.0005_dp
    type(anti_sigma_field_t) :: field_2D
    real(dp), parameter :: B_pert = 0.001_dp, M_pol_pert = 1.0_dp, &
                           N_tor_pert = 0.0_dp
    type(mock_perturbed_field_t) :: perturbed_field_2D
    type(mock_field_3d_t) :: field

    integer, parameter :: n_fieldlines = 6

    real(dp), dimension(n_fieldlines) :: xi_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    real(dp), parameter :: nfp = N_tor
    real(dp) :: iota
    type(flock_of_fieldlines_t) :: flock

    real(dp), parameter :: l_c = 1e-4, s_tor = 0.25_dp
    real(dp), parameter :: Omega_hat_zero = 0.0_dp
    real(dp), dimension(:), allocatable :: flux_modes
    type(bounce_integrals_t) :: bounce_integrals

    call field_2D%anti_sigma_field_init(M_pol, N_tor, B_0, eps_0, eps_1)
    call perturbed_field_2D%mock_perturbed_field_init(field_2D, &
                                                      M_pol_pert, &
                                                      N_tor_pert, &
                                                      B_pert)
    call field%mock_field_3d_init(perturbed_field_2D)
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

    ! Bounce integrals computed once; only Omega_hat_zero and l_c need to
    ! vary from here on.
    call setup_bounce_integrals(field, flock, s_tor, bounce_integrals)

    call compute_flux_modes(bounce_integrals, l_c, Omega_hat_zero, flux_modes)

! --------------------------------------------------------------------------
! TODO: implement assertions against flux_modes (and, if needed, rerun
! compute_flux_modes with other Omega_hat_zero/l_c values reusing the same
! bounce_integrals).
! --------------------------------------------------------------------------

end program test_precession_analytic
