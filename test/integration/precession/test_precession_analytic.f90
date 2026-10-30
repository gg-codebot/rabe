program test_precession_analytic
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use constants, only: dp, pi
    use anti_sigma_field, only: anti_sigma_field_t
    use mock_perturbed_field, only: mock_perturbed_field_t
    use mock_field_3d, only: mock_field_3d_t
    use fieldline_mod, only: flock_of_fieldlines_t
    use make_fieldline, only: make_flock_from_labels
    use precession, only: bounce_integrals_t
    use precession, only: setup_bounce_integrals
    use precession, only: compute_flux_modes
    use precession, only: get_max_valid_l_c
    use utils, only: linspace, not_same
    use myplot_module, only: myplot

    implicit none

    real(dp), parameter :: M_pol = 0.0_dp, N_tor = 2.0_dp
    real(dp), parameter :: B_0 = 1.7_dp, eps_0 = -0.0125_dp, &
                           eps_1 = -0.0005_dp
    type(anti_sigma_field_t) :: field_2D
    real(dp), parameter :: B_pert = 0.0000_dp, M_pol_pert = 1.0_dp, &
                           N_tor_pert = 0.0_dp
    type(mock_perturbed_field_t) :: perturbed_field_2D
    type(mock_field_3d_t) :: field

    integer, parameter :: n_fieldlines = 3

    real(dp), dimension(n_fieldlines) :: xi_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    real(dp), parameter :: nfp = N_tor
    real(dp) :: iota
    type(flock_of_fieldlines_t) :: flock

    real(dp), parameter :: l_c = 5.432_dp, s_tor = 0.25_dp
    integer, parameter :: n_case = 100
    integer :: case
    real(dp), dimension(n_case) :: Omega_hat
    real(dp), dimension(n_case) :: flux_1_analytic, flux_1
    real(dp), dimension(n_case) :: A_1
    real(dp) :: J_1
    real(dp), dimension(:), allocatable :: flux_modes
    type(bounce_integrals_t) :: bounce_integrals

    integer :: idx, n_eta
    real(dp), dimension(:), allocatable :: bounce_time
    real(dp), dimension(:), allocatable :: barely_trapped_tau
    real(dp), dimension(:), allocatable :: deeply_trapped_tau
    real(dp) :: shift
    logical, dimension(:), allocatable :: is_barely_trapped
    logical, dimension(:), allocatable :: is_deeply_trapped
    ! eta_tol: fraction of the trapped eta range at each end in which the
    ! deeply and barely trapped limits are compared (only valid near the ends)
    ! further in (larger eta_tol) approximations break down -> need bigger retol
    real(dp), parameter :: retol = 2e-3_dp, eta_tol = 1e-2_dp
    real(dp) :: eta_c, eta_t
    ! the asymptotic flux is valid for A_m >> 1, i.e. once the boundary layer
    ! at eta_c is thin compared to the trapped region; next order corrections
    ! are O(1/ln(A_m)), which sets the tolerance
    real(dp), parameter :: A_min = 200.0_dp
    real(dp) :: flux_tol
    ! above Omega_max the boundary layer at eta_c is not resolved by the grid
    real(dp) :: Omega_max

    logical :: test_failed

    logical :: should_plot = .false.
    type(myplot) :: plt, plt2, plt3

    test_failed = .false.

    call field_2D%anti_sigma_field_init(M_pol, N_tor, B_0, eps_0, eps_1)
    call perturbed_field_2D%mock_perturbed_field_init(field_2D, &
                                                      M_pol_pert, &
                                                      N_tor_pert, &
                                                      B_pert)
    call field%mock_field_3d_init(perturbed_field_2D)
    ! analytic J_1 assumes theta = const along a well
    field%iota = 0.0_dp
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

    ! Bounce integrals computed once
    call setup_bounce_integrals(field, flock, s_tor, bounce_integrals)

    call linspace(1.0_dp, real(n_case, kind=dp)/3.0_dp, n_case, Omega_hat)
    Omega_hat = 2.0_dp**Omega_hat
    J_1 = calc_J_1(eps_0, eps_1, B_0, N_tor)

    do case = 1, n_case
        print *, "Case ", case, " of ", n_case
        print *, "Omega_hat = ", Omega_hat(case)
        call compute_flux_modes(bounce_integrals, l_c, Omega_hat(case), flux_modes)
        print *, "Flux modes: ", flux_modes
        flux_1(case) = flux_modes(2)
        A_1(case) = calc_A_m(Omega_hat(case), l_c, 1.0_dp, eps_0, B_0, N_tor)
        flux_1_analytic(case) = calc_flux_1_analytic(A_1(case), J_1, eps_0, B_0)
    end do

    Omega_max = get_max_valid_l_c(bounce_integrals, 1.0_dp)/l_c
    print *, "Omega_max = ", Omega_max
    do case = 1, n_case
        if (A_1(case) < A_min) cycle
        if (Omega_hat(case) > Omega_max) cycle
        flux_tol = 1.0_dp/log(A_1(case))
        if (not_same(flux_1(case), flux_1_analytic(case), &
                     reltol_in=flux_tol, abstol_in=0.0_dp)) then
            print *, "-------------------------------------------------------------"
            print *, "test_precession_analytic failed: flux mode 1"
            print *, "Omega_hat = ", Omega_hat(case), "A_1 = ", A_1(case)
            print *, "found: ", flux_1(case)
            print *, "expected: ", flux_1_analytic(case)
            print *, "ratio: ", flux_1(case)/flux_1_analytic(case)
            test_failed = .true.
        end if
    end do

    n_eta = size(bounce_integrals%eta)
    allocate (bounce_time(n_eta), is_barely_trapped(n_eta), is_deeply_trapped(n_eta))
    eta_c = 1.0_dp/(B_0*(1.0_dp + abs(eps_0)))
    eta_t = 1.0_dp/(B_0*(1.0_dp - abs(eps_0) + abs(eps_1)))
    is_barely_trapped = abs(bounce_integrals%eta - eta_c) < eta_tol*(eta_t - eta_c)
    is_deeply_trapped = abs(bounce_integrals%eta - eta_t) < eta_tol*(eta_t - eta_c)

    bounce_time = bounce_integrals%avg_bounce_time_weighted/(2.0_dp*bounce_integrals%t)
    do idx = 1, size(bounce_time)
        if (bounce_integrals%t(idx) < 1e-14) then
            bounce_time(idx) = 0.0_dp
            is_barely_trapped(idx) = .false.
            is_deeply_trapped(idx) = .false.
        end if
    end do
   barely_trapped_tau = calc_barely_trapped_tau(bounce_integrals%eta, eps_0, B_0, N_tor)
   deeply_trapped_tau = calc_deeply_trapped_tau(bounce_integrals%eta, eps_0, B_0, N_tor)
    ! analytical estimate it given up to a constant shift
    idx = n_eta - 1
    shift = bounce_time(idx) - barely_trapped_tau(idx)
    barely_trapped_tau = barely_trapped_tau + shift

    if (not_same(pack(bounce_time, is_deeply_trapped), &
                 pack(deeply_trapped_tau, is_deeply_trapped), &
                 reltol_in=retol, abstol_in=0.0_dp)) then
        print *, "-------------------------------------------------------------"
        print *, "test_precession_analytic failed: deeply trapped bounce time"
        print *, "found: ", pack(bounce_time, is_deeply_trapped)
        print *, "expected: ", pack(deeply_trapped_tau, is_deeply_trapped)
        print *, "max relative error: ", maxval(abs(1.0_dp - &
                                                  pack(bounce_time, is_deeply_trapped) &
                                          /pack(deeply_trapped_tau, is_deeply_trapped)))
        test_failed = .true.
    end if
    if (not_same(pack(bounce_time, is_barely_trapped), &
                 pack(barely_trapped_tau, is_barely_trapped), &
                 reltol_in=retol, abstol_in=0.0_dp)) then
        print *, "-------------------------------------------------------------"
        print *, "test_precession_analytic failed: barely trapped bounce time"
        print *, "found: ", pack(bounce_time, is_barely_trapped)
        print *, "expected: ", pack(barely_trapped_tau, is_barely_trapped)
        print *, "max relative error: ", maxval(abs(1.0_dp - &
                                                  pack(bounce_time, is_barely_trapped) &
                                          /pack(barely_trapped_tau, is_barely_trapped)))
        test_failed = .true.
    end if

    ! I_j is not available at eta_c itself; at the next grid point the offset
    ! is O(x ln x) with x = eta - eta_c, well below retol
    idx = n_eta - 1
    if (not_same(bounce_integrals%avg_I_j(idx), calc_I_B(eps_0, B_0, N_tor), &
                 reltol_in=retol, abstol_in=0.0_dp)) then
        print *, "-------------------------------------------------------------"
        print *, "test_precession_analytic failed: I_B"
        print *, "found: ", bounce_integrals%avg_I_j(idx)
        print *, "expected: ", calc_I_B(eps_0, B_0, N_tor)
        test_failed = .true.
    end if

    if (should_plot) then
        call plot_flux(plt, plt2, Omega_hat, flux_1, flux_1_analytic)
        call plt%show()
        call plt2%show()

        call plt3%initialize(xlabel="$\sqrt{\eta - \eta_c}$", &
                             ylabel="bounce time", &
                             legend=.true., &
                             title="bounce time")
        call plt3%add_plot(bounce_integrals%t, &
                           bounce_time, &
                           label="numeric", linestyle="o-", linewidth=2)
        call plt3%add_plot(bounce_integrals%t, &
                           deeply_trapped_tau, &
                           label="deeply trapped", linestyle="-", linewidth=2)
        call plt3%add_plot(bounce_integrals%t, &
                           barely_trapped_tau, &
                           label="barely trapped", linestyle="-", linewidth=2)
        call plt3%show()
    end if

    deallocate (bounce_time)

    if (test_failed) error stop

contains

    subroutine plot_flux(plt, plt_ratio, Omega_hat, flux, flux_analytic)
        type(myplot), intent(inout) :: plt, plt_ratio
        real(dp), dimension(:), intent(in) :: Omega_hat, flux, flux_analytic

        call plt%initialize(xlabel="$\hat\Omega$", &
                            ylabel="flux mode 1", &
                            legend=.true., &
                            title="flux mode 1 vs $\hat\Omega$")
        call plt%add_plot(Omega_hat, flux, label="numeric", &
                          xscale="log", yscale="log", linestyle="o", markersize=4)
        call plt%add_plot(Omega_hat, flux_analytic, label="analytic", &
                          xscale="log", yscale="log", linestyle="-", linewidth=2)
        call plt_ratio%initialize(xlabel="$\hat\Omega$", &
                                  ylabel="numeric / analytic", &
                                  legend=.true., &
                                  title="flux mode 1 ratio")
        call plt_ratio%add_plot(Omega_hat, flux/flux_analytic, label="ratio", &
                                xscale="log", linestyle="o-", linewidth=2)
    end subroutine plot_flux

    function calc_A_m(Omega_hat, l_c, m, eps_0, B_0, N_tor)
        real(dp), intent(in) :: Omega_hat, l_c, m, eps_0, B_0, N_tor
        real(dp) :: calc_A_m
        real(dp) :: eta_b, d2B_dphi2, I_B

        eta_b = abs(1.0_dp/((1.0_dp + abs(eps_0))*B_0))
        I_B = calc_I_B(eps_0, B_0, N_tor)
        d2B_dphi2 = abs(B_0*eps_0*N_tor**2.0_dp)
        calc_A_m = sqrt(2.0_dp*m*l_c*Omega_hat*eta_b**2.0_dp/I_B)* &
                   (2.0_dp*eta_b*d2B_dphi2)**(-0.25_dp)
        if (ieee_is_nan(calc_A_m)) then
            print *, "error: A_m is NaN."
            error stop
        end if
        if (calc_A_m < 1.0_dp) then
            print *, "error: approximated formula is only valid for A_m > 1.0"
            print *, "A_m = ", calc_A_m
            print *, "Omega_hat = ", Omega_hat
            print *, "l_c = ", l_c
            print *, "m = ", m
            error stop
        end if
    end function calc_A_m

    !> field line average of I_j at eta_c (independent of iota to O(eps_1^2))
    function calc_I_B(eps_0, B_0, N_tor) result(I_B)
        real(dp), intent(in) :: eps_0, B_0, N_tor
        real(dp) :: I_B

        I_B = 4.0_dp*sqrt(2.0_dp*abs(eps_0)/(1.0_dp + abs(eps_0)))/N_tor/B_0**2.0_dp
        I_B = I_B*(1.0_dp + 2.0_dp/3.0_dp*abs(eps_0))
    end function calc_I_B

    !> m=1 mode at eta_c, J_1 = -2/3 d_eta I_V,1 (cf. test_radial_drift),
    !! valid for iota = 0
    function calc_J_1(eps_0, eps_1, B_0, N_tor)
        real(dp), intent(in) :: eps_0, eps_1, B_0, N_tor
        real(dp) :: calc_J_1

        calc_J_1 = 2.0_dp*sqrt(2.0_dp*abs(eps_0)/(1.0_dp + abs(eps_0)))/N_tor/B_0* &
                   eps_1/eps_0*(1.0_dp + 3.0_dp*abs(eps_0))
        if (ieee_is_nan(calc_J_1)) then
            print *, "error: calc_J_1 is NaN."
            error stop
        end if
    end function calc_J_1

    function calc_flux_1_analytic(A_1, J_1, eps_0, B_0) result(flux_1)
        real(dp), intent(in) :: A_1, J_1, eps_0, B_0
        real(dp) :: flux_1
        real(dp) :: eta_b

        eta_b = abs(1.0_dp/((1.0_dp + abs(eps_0))*B_0))
        flux_1 = sqrt(2.0_dp)*eta_b*J_1/(A_1*sqrt(abs(log(A_1))))
        if (ieee_is_nan(flux_1)) then
            print *, "error: flux_1 is NaN."
            error stop
        end if
    end function calc_flux_1_analytic

    function calc_deeply_trapped_tau(eta, eps_0, B_0, N_tor) result(tau)
        real(dp), dimension(:), intent(in) :: eta
        real(dp), intent(in) :: eps_0, B_0, N_tor

        real(dp), dimension(size(eta)) :: tau
        real(dp), dimension(size(eta)) :: k_squared
        real(dp) :: d2B_dphi2, B_min

        d2B_dphi2 = abs(B_0*eps_0*N_tor**2.0_dp)
        B_min = B_0*(1.0_dp - abs(eps_0))
        k_squared = (1.0_dp - B_min*eta)/(2.0_dp*eta*abs(B_0*eps_0))
        ! deeply trapped particles see B ~ B_min
        tau = 4.0_dp*pi*sqrt(0.5_dp/(d2B_dphi2*eta))*(1.0_dp + k_squared/4.0_dp) &
              /B_min

    end function calc_deeply_trapped_tau

    function calc_barely_trapped_tau(eta, eps_0, B_0, N_tor) result(tau)
        real(dp), dimension(:), intent(in) :: eta
        real(dp), intent(in) :: eps_0, B_0, N_tor

        real(dp), dimension(size(eta)) :: tau
        real(dp), dimension(size(eta)) :: k_squared
        real(dp) :: d2B_dphi2, eta_c

        d2B_dphi2 = abs(B_0*eps_0*N_tor**2.0_dp)
        eta_c = 1.0_dp/(B_0*(1.0_dp + abs(eps_0)))
        ! barely trapped particles see B ~ B_max = 1/eta_c
        tau = -4.0_dp*sqrt(0.5_dp/(d2B_dphi2*eta))*log((eta - eta_c)/eta_c) &
              *eta_c

    end function calc_barely_trapped_tau

end program test_precession_analytic
