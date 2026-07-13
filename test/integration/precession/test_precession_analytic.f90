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
    use utils, only: linspace
    use myplot_module, only: myplot

    implicit none

    real(dp), parameter :: M_pol = 0.0_dp, N_tor = 1.0_dp
    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = -0.0125_dp, &
                           eps_1 = -0.0005_dp
    type(anti_sigma_field_t) :: field_2D
    real(dp), parameter :: B_pert = 0.0001_dp, M_pol_pert = 1.0_dp, &
                           N_tor_pert = 0.0_dp
    type(mock_perturbed_field_t) :: perturbed_field_2D
    type(mock_field_3d_t) :: field

    integer, parameter :: n_fieldlines = 6

    real(dp), dimension(n_fieldlines) :: xi_0
    real(dp), dimension(n_fieldlines + 1) :: temp
    real(dp), parameter :: nfp = N_tor
    real(dp) :: iota
    type(flock_of_fieldlines_t) :: flock

    real(dp), parameter :: l_c = 1.0_dp, s_tor = 0.25_dp
    integer, parameter :: n_case = 10
    integer :: case
    real(dp), dimension(n_case) :: Omega_hat
    real(dp), dimension(n_case) :: flux_1_analytic, flux_1
    real(dp) :: A_1, J_1
    real(dp), dimension(:), allocatable :: flux_modes
    type(bounce_integrals_t) :: bounce_integrals

    logical :: should_plot = .true.
    type(myplot) :: plt

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

    call linspace(1.0_dp, real(n_case, kind=dp), n_case, Omega_hat)
    Omega_hat = 2.0_dp**Omega_hat
    J_1 = calc_J_1(eps_0, eps_1, B_0, N_tor)

    do case = 1, n_case
        print *, "Case ", case, " of ", n_case
        print *, "Omega_hat = ", Omega_hat(case)
        call compute_flux_modes(bounce_integrals, l_c, Omega_hat(case), flux_modes)
        print *, "Flux modes: ", flux_modes
        flux_1(case) = flux_modes(2)
        A_1 = calc_A_m(Omega_hat(case), l_c, 1.0_dp, eps_0, B_0, N_tor)
        flux_1_analytic(case) = calc_flux_1_analytic(A_1, J_1, eps_0, B_0)
    end do

    if (should_plot) then
        call plt%initialize(xlabel="log(Omega_hat)", &
                            ylabel="Flux mode 1", &
                            legend=.true., &
                            title="Flux mode 1 vs Omega_hat")
        call plt%add_plot(log10(Omega_hat), flux_1, label="numerical", &
                          xscale="log", linestyle="-", linewidth=2)
        call plt%add_plot(log10(Omega_hat), flux_1_analytic, label="analytic", &
                          xscale="log", linestyle="--", linewidth=2)
        call plt%show()
    end if

contains

    function calc_A_m(Omega_hat, l_c, m, eps_0, B_0, N_tor)
        real(dp), intent(in) :: Omega_hat, l_c, m, eps_0, B_0, N_tor
        real(dp) :: calc_A_m
        real(dp) :: eta_b, d2B_dphi2, I_B

        eta_b = abs(1.0_dp/((1.0_dp + abs(eps_0))*B_0))
        I_B = 4.0_dp*sqrt(2.0_dp*abs(eps_0)/(1.0_dp + abs(eps_0)))/N_tor/B_0**2.0_dp
        I_B = I_B*(1.0_dp + 2.0_dp/3.0_dp*abs(eps_0))
        d2B_dphi2 = abs(B_0*eps_0*N_tor**2.0_dp)
        calc_A_m = sqrt(pi*m*l_c*Omega_hat*eta_b**2.0_dp/I_B)* &
                   (2.0_dp*eta_b*d2B_dphi2)**0.25_dp
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

    function calc_J_1(eps_0, eps_1, B_0, N_tor)
        real(dp), intent(in) :: eps_0, eps_1, B_0, N_tor
        real(dp) :: calc_J_1
        real(dp) :: eta_b, I_B

        I_B = 4.0_dp*sqrt(2.0_dp*abs(eps_0)/(1.0_dp + abs(eps_0)))/N_tor/B_0**2.0_dp
        I_B = I_B*(1.0_dp + 2.0_dp/3.0_dp*abs(eps_0))
        calc_J_1 = -0.25_dp*eps_1/eps_0*(1.0_dp + 16.0_dp/3.0_dp*abs(eps_0))
        calc_J_1 = calc_J_1*I_B*B_0
        if (ieee_is_nan(calc_J_1)) then
            print *, "error: calc_J_1 is NaN."
            error stop
        end if
    end function calc_J_1

    function calc_flux_1_analytic(A_1, J_1, eps_0, B_0) result(flux_1)
        real(dp), intent(in) :: A_1, J_1, eps_0, B_0
        real(dp) :: flux_1
        real(dp) :: eta_b, I_B

        eta_b = abs(1.0_dp/((1.0_dp + abs(eps_0))*B_0))
        flux_1 = sqrt(2.0_dp)*eta_b*J_1/(A_1*sqrt(abs(log(A_1))))
        if (ieee_is_nan(flux_1)) then
            print *, "error: flux_1 is NaN."
            error stop
        end if
    end function calc_flux_1_analytic

end program test_precession_analytic
