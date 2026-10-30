program test_zero_precession
    use constants, only: dp, pi
    use anti_sigma_field, only: anti_sigma_field_t
    use mock_perturbed_field, only: mock_perturbed_field_t
    use mock_field_3d, only: mock_field_3d_t
    use fieldline_mod, only: flock_of_fieldlines_t
    use make_fieldline, only: make_flock_from_labels
    use precession, only: bounce_integrals_t
    use precession, only: setup_bounce_integrals
    use precession, only: compute_flux_modes
    use precession_reference, only: integrate_radial_drift
    use grid_mod, only: fieldline_for_precession_t
    use grid_mod, only: set_splines
    use fourier, only: real_ft
    use utils, only: not_same, linspace

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
    real(dp), parameter :: odetol = 1e-8
    real(dp) :: tol
    type(bounce_integrals_t) :: bounce_integrals
    real(dp), dimension(:), allocatable :: flux_modes
    integer :: n_modes, idx
    logical :: test_failed

    real(dp) :: expected, found

    ! numerical check variables
    type(fieldline_for_precession_t), dimension(:), allocatable :: fieldlines
    real(dp), dimension(n_fieldlines) :: radial_drift_integrals
    integer, parameter :: n_modes_ft = n_fieldlines/2 + 1
    real(dp), dimension(n_modes_ft) :: numerical_cos, numerical_sin

    test_failed = .false.

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

    call setup_bounce_integrals(field, flock, s_tor, bounce_integrals, &
                                fieldlines_out=fieldlines)
    call compute_flux_modes(bounce_integrals, l_c, Omega_hat_zero, flux_modes)

    n_modes = n_fieldlines/2 + 1

    ! reference: same bounce samples, integrated exactly per fieldline
    do idx = 1, n_fieldlines
        call set_splines(fieldlines(idx)%grid)
        call integrate_radial_drift(fieldlines(idx)%grid, &
                                    fieldlines(idx)%fieldline_t, &
                                    radial_drift_integrals(idx))
    end do

    call real_ft(xi_0, radial_drift_integrals, &
                 numerical_cos, numerical_sin)

    ! assuming odeint completes each integration with
    ! a single step i.e. lower limit of worst case error
    ! -> a failing test might be due to too strict tolerances
    tol = odetol*fieldlines(1)%grid%n
    do idx = 2, n_modes
        expected = numerical_sin(idx)
        found = flux_modes(idx)
        if (not_same(found, expected, reltol_in=tol, &
                     abstol_in=0.0_dp)) then
            print *, "---------------------------------------------------"
            print *, "test_zero_precession failed: mode", idx
            print *, "found: ", found
            print *, "expected: ", expected
            print *, "relative error: ", abs(found - expected)/abs(expected)
            test_failed = .true.
        end if
    end do

    deallocate (fieldlines)

    if (test_failed) error stop

end program test_zero_precession
