program test_radial_drift
    use constants, only: dp, pi
    use anti_sigma_field, only: anti_sigma_field_t
    use mock_perturbed_field, only: mock_perturbed_field_t
    use mock_field_3d, only: mock_field_3d_t
    use fieldline_mod, only: flock_of_fieldlines_t
    use make_fieldline, only: make_flock_from_labels
    use precession, only: bounce_integrals_t
    use precession, only: setup_bounce_integrals
    use utils, only: linspace, not_same
    use myplot_module, only: myplot
    use grid_mod, only: fieldline_for_precession_t

    implicit none

    real(dp), parameter :: M_pol = 0.0_dp, N_tor = 2.0_dp
    real(dp), parameter :: B_0 = 2.0_dp, eps_0 = -0.0125_dp, &
                           eps_1 = -0.0001_dp
    type(anti_sigma_field_t) :: field_2D
    real(dp), parameter :: B_pert = 0.0000_dp, M_pol_pert = 1.0_dp, &
                           N_tor_pert = 0.0_dp
    type(mock_perturbed_field_t) :: perturbed_field_2D
    type(mock_field_3d_t) :: field

    integer, parameter :: n_fieldlines = 80

    real(dp), dimension(n_fieldlines) :: xi_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    real(dp), parameter :: nfp = N_tor
    real(dp) :: iota
    type(flock_of_fieldlines_t) :: flock

    real(dp), parameter :: s_tor = 0.25_dp
    type(bounce_integrals_t) :: bounce_integrals

    integer :: idx
    integer, parameter :: n_eta = 20
    integer :: idx_eta
    real(dp), dimension(:), allocatable :: analytic_drift, drift
    real(dp) :: drift_1_analytic, drift_1
    real(dp) :: dAtheta_dr, B_theta_covar, B_phi_covar
    real(dp), dimension(:, :), allocatable :: modes_I_B, modes_dIV_deta
    real(dp) :: reltol, abstol, max_error, delta_eta, delta_trapped
    real(dp) :: B_min, B_varphi, B_vartheta
    type(fieldline_for_precession_t), dimension(:), allocatable :: fieldlines

    logical :: test_failed

    logical :: should_plot = .true.
    type(myplot) :: plt
    integer :: n_plot_points = 100
    real(dp), dimension(:), allocatable :: xi_0_plot
    real(dp), dimension(:), allocatable :: I_j, dIV_deta, dIj_dxi0, int_dIV_deta
    real(dp), dimension(n_eta) :: sum_modes
    character(len=256) :: label

    test_failed = .false.

    call field_2D%anti_sigma_field_init(M_pol, N_tor, B_0, eps_0, eps_1)
    call perturbed_field_2D%mock_perturbed_field_init(field_2D, &
                                                      M_pol_pert, &
                                                      N_tor_pert, &
                                                      B_pert)
    call field%mock_field_3d_init(perturbed_field_2D)
    field%iota = 0.01_dp
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

    call setup_bounce_integrals(field, flock, s_tor, bounce_integrals, n_eta=n_eta, fieldlines_out=fieldlines)

    call field%get_covariant_components(B_theta_covar, B_phi_covar)
    ! Mock field is a right-handed system (sqrtg>0), so dAtheta_dr = dpsi_dr
    ! and radial coordinate is s_tor so that dpsi_dr = psi_tor_edge
    dAtheta_dr = field%psi_tor_edge
    drift_1_analytic = barely_trapped_drift_mode(eps_0, eps_1, B_0, dAtheta_dr, &
                                                B_theta_covar, B_phi_covar, iota, N_tor)
    drift_1 = bounce_integrals%radial_drift_sin(2,n_eta-1)/(2.0_dp*bounce_integrals%t(n_eta-1))
    ! numerical drift is normalized to B*rho_L i.e. units of 1/T
    ! analytic drift is normalized to rho_L*B/B_0 i.e. no units, for comparison
    drift_1 = drift_1*B_0

    ! Numerical drift not available exactly at trapped-passing boundary, so
    ! error is poluted by difference in eta of those points
    ! -> can be further tighted by higher n_eta
    delta_eta = abs(bounce_integrals%eta(n_eta) - bounce_integrals%eta(n_eta - 1))
    delta_trapped = abs(bounce_integrals%eta(n_eta) - bounce_integrals%eta(1))
    reltol = eps_0**2 + delta_eta/delta_trapped*log(delta_trapped/delta_eta)
    print *, "reltol = ", reltol
    if (not_same(drift_1, drift_1_analytic, &
                 reltol_in=reltol, abstol_in=0.0_dp)) then
        print *, "test_precession_analytic failed: radial drift m=1"
        print *, "found: ", drift_1
        print *, "expected: ", drift_1_analytic
        print *, "relative error: ", abs(drift_1 - drift_1_analytic)/ &
            abs(drift_1_analytic)
        test_failed = .true.
    end if

    do idx = 1, n_eta
        sum_modes(idx) = sum(abs(bounce_integrals%radial_drift_sin(:, idx)))
        sum_modes(idx) = sum_modes(idx) - &
                         abs(bounce_integrals%radial_drift_sin(2, idx))
        sum_modes(idx) = sum_modes(idx)/(2.0_dp*bounce_integrals%t(idx))*B_0
    end do
    sum_modes = sum_modes/abs(drift_1)
    if (not_same(sum_modes(n_eta - 1), 0.0_dp, &
                 reltol_in=0.0_dp, abstol_in=reltol)) then
        print *, "test_precession_analytic failed: sum radial drift m>0 at eta_c"
        print *, "found: ", sum_modes
        print *, "expected: ", 0.0_dp
        print *, "absolute error: ", sum_modes(n_eta - 1)
        test_failed = .true.
    end if

    allocate (int_dIV_deta(n_fieldlines))
    allocate (dIV_deta(n_fieldlines))
    allocate (I_j(n_fieldlines))
    allocate (dIj_dxi0(n_fieldlines))
    do idx = 1, n_eta - 1
        call get_dIV_deta_and_antiderivative(fieldlines, idx, B_0, dAtheta_dr, &
                               B_theta_covar, B_phi_covar, iota, dIV_deta, int_dIV_deta)
        call get_I_j_and_derivative(fieldlines, idx, xi_0, I_j, dIj_dxi0)
        abstol = maxval(abs(dIV_deta))*reltol
        if (not_same(dIV_deta, 1.5_dp*B_0*dIj_dxi0, &
                     reltol_in=0.0_dp, abstol_in=abstol)) then
            print *, "test_precession_analytic failed: dIV_deta vs dIj_dxi0"
            print *, "eta_idx = ", idx, "eta = ", bounce_integrals%eta(idx)
            max_error = maxval(abs(dIV_deta - 1.5_dp*B_0*dIj_dxi0))
            print *, "max rel error: ", max_error/maxval(abs(dIV_deta))
            test_failed = .true.
        end if
        abstol = maxval(abs(int_dIV_deta))*reltol
        if (not_same(int_dIV_deta, &
                     1.5_dp*B_0*(I_j - I_j(1)), &
                     reltol_in=0.0_dp, abstol_in=abstol)) then
            print *, "test_precession_analytic failed: int_dIV_deta vs I_j"
            print *, "eta_idx = ", idx, "eta = ", bounce_integrals%eta(idx)
            max_error = maxval(abs(int_dIV_deta - 1.5_dp*B_0*(I_j - I_j(1))))
            print *, "max abs error: ", max_error
            print *, "max abs int_dIV_deta: ", maxval(abs(int_dIV_deta))
            print *, "max rel error: ", max_error/maxval(abs(int_dIV_deta))
            test_failed = .true.
        end if
    end do

    if (should_plot) then
        ! plot barely trapped radial drift
        call plt%initialize(xlabel="$\xi_0$", &
                            ylabel="drift [1]", &
                            legend=.true., &
                            title="normalized radial drift for barely trapped", &
                            fontsize=20)
        allocate (drift(n_plot_points))
        allocate (xi_0_plot(n_plot_points))
        call linspace(0.0_dp, 2.0_dp*pi, n_plot_points, xi_0_plot)
        do idx = n_eta - 9, n_eta - 1
            write (label, '(A,F6.3)') "numeric mode, eta = ", bounce_integrals%eta(idx)

            drift = bounce_integrals%radial_drift_sin(2, idx)*sin(xi_0_plot)
            drift = drift/(2.0_dp*bounce_integrals%t(idx))*B_0
            call plt%add_plot(xi_0_plot, &
                              drift, &
                              label=label, linestyle="--", linewidth=2)
        end do
        allocate (analytic_drift(n_plot_points))
        analytic_drift = drift_1_analytic*sin(xi_0_plot)
        call plt%add_plot(xi_0_plot, &
                          analytic_drift, &
                          label="analytic", linestyle="-", linewidth=4)
        deallocate (drift)
        allocate (drift(n_fieldlines))
        do idx = 1, n_fieldlines
            drift(idx) = fieldlines(idx)%grid%radial_drift_weighted(n_eta - 1)
            drift(idx) = drift(idx)/(2.0_dp*bounce_integrals%t(n_eta - 1))*B_0
        end do
        write (label, '(A,F6.3)') "numeric, eta = ", bounce_integrals%eta(n_eta - 1)
        call plt%add_plot(xi_0, &
                          drift, &
                          label=label, linestyle="s", markersize=10)

        ! plot remainder of modes m>1
        call plt%initialize(xlabel="$\sqrt{\eta - \eta_c}$", &
                            ylabel="$\sum |$drift$^m|/|$drift$^1|$ [1]", &
                            legend=.true., &
                            title="sum of modes m > 1", &
                            fontsize=20)
        call plt%add_plot(bounce_integrals%t(1:n_eta - 1), &
                          sum_modes(1:n_eta - 1), &
                          label="normalized by m=1 mode", linestyle="-", linewidth=4)

        ! plot comparison of dIV_deta to I_j
        idx_eta = n_eta - 1
        write (label, '(A,F6.3)') "comparison $\partial_\eta I_\mathcal{V}$ to $I_j$ at eta = ", &
            bounce_integrals%eta(idx_eta)
        call plt%initialize(xlabel="$\xi_0$", &
                            ylabel="drift [1]", &
                            legend=.true., &
                            title=label, &
                            fontsize=20)
        call get_dIV_deta_and_antiderivative(fieldlines, idx_eta, B_0, dAtheta_dr, &
                               B_theta_covar, B_phi_covar, iota, dIV_deta, int_dIV_deta)
        call get_I_j_and_derivative(fieldlines, idx_eta, xi_0, I_j, dIj_dxi0)
        call plt%add_plot(xi_0, &
                          1.5_dp*B_0*(I_j - I_j(1)), &
                    label="$1.5 B_0 (I_j - I_j(\xi_0=0))$", linestyle="-o", linewidth=2)
        call plt%add_plot(xi_0, &
                          dIV_deta, &
                     label="$\partial_\eta I_\mathcal{V}$", linestyle="-o", linewidth=2)
        call plt%add_plot(xi_0, &
                          1.5*B_0*dIj_dxi0, &
                    label="$1.5 B_0 \partial_{\xi_0} I_j$", linestyle="-o", linewidth=2)
        call plt%add_plot(xi_0, &
                          int_dIV_deta, &
        label="$\int d \xi_0 \partial_\eta I_\mathcal{V}$", linestyle="-o", linewidth=2)
        call plt%show()
    end if

    if (test_failed) error stop

contains

    !> compute the first mode in terms of field line label of the radial drift
    !! <Delta_r>/(B*rho_L)*B_0
    !! normalized by the reference lamor radius m*v/(q*B_0)
    function barely_trapped_drift_mode(eps_0, eps_1, B_0, dAtheta_dr, &
                                B_theta_covar, B_phi_covar, iota, N_tor) result(drift_1)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
        real(dp), intent(in) :: eps_0, eps_1, B_0, dAtheta_dr
        real(dp), intent(in) :: B_theta_covar, B_phi_covar, iota
        real(dp), intent(in) :: N_tor
        real(dp) :: drift_1

        real(dp) :: B_max, dI_V_1_deta

        B_max = abs(B_0*(1.0_dp + abs(eps_0)))
        dI_V_1_deta = -3.0_dp/B_0*sqrt(2.0_dp*abs(B_0*eps_0)/B_max)* &
                      abs(eps_1/eps_0)*(1.0_dp + 3.0_dp*abs(eps_0))/N_tor
        drift_1 = dI_V_1_deta/0.75_dp*(B_phi_covar + iota*B_theta_covar)/dAtheta_dr*B_0

        if (ieee_is_nan(drift_1)) then
            print *, "drift_1 is NaN!"
            print *, "eps_0 = ", eps_0
            print *, "eps_1 = ", eps_1
            print *, "B_0 = ", B_0
            print *, "dAtheta_dr = ", dAtheta_dr
            print *, "B_theta_covar = ", B_theta_covar
            print *, "B_phi_covar = ", B_phi_covar
            print *, "iota = ", iota
            print *, "B_max = ", B_max
            error stop
        end if

    end function barely_trapped_drift_mode

    subroutine get_dIV_deta_and_antiderivative(fieldlines, eta_idx, B_0, dAtheta_dr, &
                               B_theta_covar, B_phi_covar, iota, dIV_deta, int_dIV_deta)
        type(fieldline_for_precession_t), dimension(:), intent(in) :: fieldlines
        integer, intent(in) :: eta_idx
        real(dp), intent(in) :: B_0
        real(dp), intent(in) :: dAtheta_dr, B_theta_covar, B_phi_covar, iota
        real(dp), dimension(size(fieldlines)), intent(out) :: dIV_deta, int_dIV_deta

        integer :: idx

        real(dp), dimension(size(fieldlines)) :: drift
        do idx = 1, n_fieldlines
            drift(idx) = fieldlines(idx)%grid%normalized_radial_drift(eta_idx)*B_0
        end do
        dIV_deta = calc_dIV_deta(drift, dAtheta_dr, B_theta_covar, B_phi_covar, iota)
        int_dIV_deta = get_antiderivative(dIV_deta, xi_0)
    end subroutine get_dIV_deta_and_antiderivative

    function calc_dIV_deta(drift, dAtheta_dr, B_theta_covar, B_phi_covar, iota) result(dIV_deta)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
        real(dp), intent(in), dimension(:) :: drift
        real(dp), intent(in) :: dAtheta_dr, B_theta_covar, B_phi_covar, iota
        real(dp), dimension(size(drift)) :: dIV_deta

        dIV_deta = drift*0.75_dp/(B_phi_covar + iota*B_theta_covar)*dAtheta_dr/B_0

        if (any(ieee_is_nan(dIV_deta))) then
            print *, "dIV_deta is NaN!"
            print *, "drift = ", drift
            print *, "dAtheta_dr = ", dAtheta_dr
            print *, "B_theta_covar = ", B_theta_covar
            print *, "B_phi_covar = ", B_phi_covar
            print *, "iota = ", iota
            error stop
        end if
    end function calc_dIV_deta

    function get_antiderivative(func, xi_0) result(antiderivative)
        real(dp), intent(in), dimension(:) :: func, xi_0
        real(dp), dimension(size(func)) :: antiderivative

        integer :: idx

        antiderivative(1) = 0.0_dp
        do idx = 2, size(func)
            antiderivative(idx) = antiderivative(idx - 1) + &
                          0.5_dp*(func(idx) + func(idx - 1))*(xi_0(idx) - xi_0(idx - 1))
        end do
    end function get_antiderivative

    subroutine get_I_j_and_derivative(fieldlines, eta_idx, xi_0, I_j, dIj_dxi0)
        type(fieldline_for_precession_t), dimension(:), intent(in) :: fieldlines
        integer, intent(in) :: eta_idx
        real(dp), intent(in), dimension(:) :: xi_0
        real(dp), dimension(size(fieldlines)), intent(out) :: I_j, dIj_dxi0

        integer :: idx

        do idx = 1, size(fieldlines)
            I_j(idx) = fieldlines(idx)%grid%I_j(eta_idx)
        end do
        dIj_dxi0 = get_derivative(I_j, xi_0)
    end subroutine get_I_j_and_derivative

    function get_derivative(func, xi_0) result(derivative)
        real(dp), intent(in), dimension(:) :: func, xi_0
        real(dp), dimension(size(func)) :: derivative

        integer :: n, idx

        n = size(xi_0)
        derivative(1) = (func(2) - func(1))/(xi_0(2) - xi_0(1))
        derivative(n) = (func(n) - func(n - 1))/(xi_0(n) - xi_0(n - 1))
        do idx = 2, n - 1
       derivative(idx) = (func(idx + 1) - func(idx - 1))/(xi_0(idx + 1) - xi_0(idx - 1))
        end do
    end function get_derivative

end program test_radial_drift
