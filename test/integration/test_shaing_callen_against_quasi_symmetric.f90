program test_shaing_callen_against_quasi_symmetric
    use mock_field, only: mock_field_t
    use fieldline_mod, only: flock_of_fieldlines_t
    use make_fieldline, only: make_flock_from_labels
    use constants, only: dp, pi
    use utils, only: linspace

    use test_shaing_callen_mod, only: test_trapped_fraction_against_qs
    use test_shaing_callen_mod, only: test_calc_avg_normalized_B_squared_dphimax_dxi0
    use test_shaing_callen_mod, only: test_calc_avg_normalized_lambda_dphimax_dxi0
    use test_shaing_callen_mod, only: test_get_non_omnigenous_remainder
    use test_shaing_callen_mod, only: test_trapped_fraction_against_circular_tokamak

    use plot_quantities, only: plot_fieldlines_over_field
    use plot_quantities, only: plot_phi_max_over_xi_0

    implicit none

    integer, parameter :: n_fieldlines = 6

    real(dp), dimension(n_fieldlines) :: xi_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    type(flock_of_fieldlines_t) :: flock

    real(dp), parameter :: nfp = 10.0_dp, iota = 0.47_dp
    real(dp), parameter :: M_pol = 2.0_dp, N_tor = nfp
    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = -0.01_dp
    type(mock_field_t) :: field

    logical, parameter :: should_plot = .false.
    integer :: this
    integer, parameter :: n_eta = 10
    real(dp), dimension(n_eta) :: eta_grid

    logical :: test_failed

    test_failed = .false.

    call field%mock_field_init(M_pol, N_tor, B_0, eps_0)

    call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines + 1, temp)
    xi_0 = temp(1:n_fieldlines)

    call make_flock_from_labels(flock, &
                                xi_0, &
                                iota, &
                                field, &
                                M_pol, &
                                N_tor, &
                                nfp)

    if (should_plot) then
        call plot_fieldlines_over_field(flock%fieldlines, field)
        call plot_phi_max_over_xi_0(flock%fieldlines, flock%M_pol, flock%nfp)
    end if

    call test_trapped_fraction_against_circular_tokamak(test_failed)
    call test_trapped_fraction_against_qs(field, flock, test_failed)
    call test_calc_avg_normalized_B_squared_dphimax_dxi0(flock, test_failed)
    call test_calc_avg_normalized_lambda_dphimax_dxi0(field, flock, test_failed)
    call test_get_non_omnigenous_remainder(field, flock, test_failed)

    if (test_failed) error stop

end program test_shaing_callen_against_quasi_symmetric
