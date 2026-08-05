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

    implicit none

    real(dp), parameter :: M_pol = 0.0_dp, N_tor = 2.0_dp
    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = -0.0225_dp, &
                           eps_1 = -0.0001_dp
    type(anti_sigma_field_t) :: field_2D
    real(dp), parameter :: B_pert = 0.0000_dp, M_pol_pert = 1.0_dp, &
                           N_tor_pert = 0.0_dp
    type(mock_perturbed_field_t) :: perturbed_field_2D
    type(mock_field_3d_t) :: field

    integer, parameter :: n_fieldlines = 10

    real(dp), dimension(n_fieldlines) :: xi_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    real(dp), parameter :: nfp = N_tor
    real(dp) :: iota
    type(flock_of_fieldlines_t) :: flock

    real(dp), parameter :: s_tor = 0.25_dp
    type(bounce_integrals_t) :: bounce_integrals

    integer :: idx
    integer, parameter :: n_eta = 20
    real(dp), dimension(:), allocatable :: analytic_drift, drift
    real(dp) :: drift_1_analytic, drift_1
    real(dp) :: dAtheta_dr, B_theta_covar, B_phi_covar
    real(dp), dimension(:, :), allocatable :: modes_I_B, modes_dIV_deta
    real(dp) :: reltol, delta_eta, delta_trapped
    real(dp) :: B_min, B_varphi, B_vartheta

    logical :: test_failed

    logical :: should_plot = .true.
    type(myplot) :: plt
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

    call setup_bounce_integrals(field, flock, s_tor, bounce_integrals, n_eta=n_eta)

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
                 reltol_in=-reltol, abstol_in=0.0_dp)) then
        print *, "test_precession_analytic failed: radial drift"
        print *, "found: ", drift_1
        print *, "expected: ", drift_1_analytic
      print *, "relative error: ", abs(drift_1 - drift_1_analytic)/abs(drift_1_analytic)
        test_failed = .true.
    end if

    if (should_plot) then
        call plt%initialize(xlabel="$\xi_0$", &
                            ylabel="drift [1]", &
                            legend=.true., &
                            title="normalized radial drift for barely trapped", &
                            fontsize=20)
        allocate (drift(n_fieldlines))
        do idx = n_eta - 9, n_eta - 1
            write (label, '(A,F6.3)') "numeric, eta = ", bounce_integrals%eta(idx)

            drift = bounce_integrals%radial_drift_sin(2, idx)*sin(xi_0)
            drift = drift/(2.0_dp*bounce_integrals%t(idx))
            call plt%add_plot(xi_0, &
                              drift, &
                              label=label, linestyle="-o", linewidth=2)
        end do

        allocate (analytic_drift(n_fieldlines))
        analytic_drift = drift_1_analytic*sin(xi_0)
        call plt%add_plot(xi_0, &
                          analytic_drift, &
                          label="analytic", linestyle="-", linewidth=4)
        call plt%show()
        call plt%initialize(xlabel="$\xi_0$", &
                            ylabel="drift [1]", &
                            legend=.true., &
                            title="other modes", &
                            fontsize=20)
        write (label, '(A,F6.3)') "numeric, eta = ", bounce_integrals%eta(idx)

        drift = sum(abs(bounce_integrals%radial_drift_sin), dim=2)
        drift = drift/(2.0_dp*bounce_integrals%t(idx))
        call plt%add_plot(xi_0, &
                          drift, &
                          label=label, linestyle="-o", linewidth=2)
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

    function calc_dIV_deta(drift, dAtheta_dr, B_theta_covar, B_phi_covar, iota) result(dIV_deta)
        real(dp), intent(in), dimension(:) :: drift
        real(dp), intent(in) :: dAtheta_dr, B_theta_covar, B_phi_covar, iota
        real(dp), dimension(size(drift)) :: dIV_deta

        dIV_deta = drift*0.75_dp/(B_phi_covar + iota*B_theta_covar)*dAtheta_dr/B_0

    end function calc_dIV_deta

end program test_radial_drift
