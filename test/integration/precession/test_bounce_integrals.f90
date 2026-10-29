program test_precession_analytic
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

    real(dp), parameter :: M_pol = 0.0_dp, N_tor = 3.0_dp
    real(dp), parameter :: B_0 = 1.5_dp, eps_0 = -0.0125_dp, &
                           eps_1 = -0.0000_dp
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

    real(dp), parameter :: s_tor = 0.25_dp
    type(bounce_integrals_t) :: bounce_integrals

    integer :: idx, n_eta
    real(dp), dimension(:), allocatable :: tau, analytic_tau
    real(dp), dimension(:), allocatable :: deeply_trapped_tau
    logical, dimension(:), allocatable :: is_valid
    real(dp), parameter :: retol = abs(eps_0)
    real(dp) :: B_min, B_varphi, B_vartheta

    logical :: test_failed

    logical :: should_plot = .false.
    type(myplot) :: plt

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

    call setup_bounce_integrals(field, flock, s_tor, bounce_integrals)

    n_eta = size(bounce_integrals%eta)
    allocate (tau(n_eta))
    allocate (analytic_tau(n_eta))
    allocate (deeply_trapped_tau(n_eta))
    allocate (is_valid(n_eta))
    is_valid = .true.

    ! tau =  bounce_time * v / (B_varphi + iota*B_vartheta) i.e. [1/T]
    tau = bounce_integrals%avg_bounce_time_weighted/(2.0_dp*bounce_integrals%t)
    call field%get_covariant_components(B_vartheta, B_varphi)
    tau = tau/(B_varphi + iota*B_vartheta)
    do idx = 1, size(tau)
        if (bounce_integrals%t(idx) < 1e-14) then
            tau(idx) = 0.0_dp
            is_valid(idx) = .false.
        end if
    end do
    B_min = 1.0_dp/maxval(bounce_integrals%eta)
    analytic_tau = calc_tau(bounce_integrals%eta, eps_0, B_0, B_min, N_tor)
   deeply_trapped_tau = calc_deeply_trapped_tau(bounce_integrals%eta, eps_0, B_0, N_tor)

    if (not_same(pack(tau, is_valid), pack(analytic_tau, is_valid), &
                 reltol_in=retol, abstol_in=0.0_dp)) then
        print *, "test_precession_analytic failed: bounce time"
        print *, "found: ", pack(tau, is_valid)
        print *, "expected: ", pack(analytic_tau, is_valid)
        print *, "max relative error: ", maxval(abs(pack(tau, is_valid) - &
                                                    pack(analytic_tau, is_valid))/ &
                                                abs(pack(tau, is_valid)))
        test_failed = .true.
    end if

    if (should_plot) then
        call plt%initialize(xlabel="$\sqrt{\eta - \eta_c}$", &
                            ylabel="$\tau$ [1/T]", &
                            legend=.true., &
                            title="bounce time")
        call plt%add_plot(bounce_integrals%t, &
                          tau, &
                          label="numeric", linestyle="o-", linewidth=2)
        call plt%add_plot(bounce_integrals%t, &
                          analytic_tau, &
                          label="analytic", linestyle="-", linewidth=2)
        call plt%add_plot(bounce_integrals%t, &
                          deeply_trapped_tau, &
                          label="deeply trapped approx", linestyle="--", linewidth=2)
        call plt%show()
    end if

    deallocate (tau)

    if (test_failed) error stop

contains

    function calc_tau(eta, eps_0, B_0, B_min, N_tor) result(tau)
        use constants, only: machine_eps
        real(dp), intent(in), dimension(:) :: eta
        real(dp), intent(in) :: eps_0, B_0, B_min, N_tor

        real(dp), dimension(size(eta)) :: tau
        real(dp), dimension(size(eta)) :: k_squared

        integer :: idx
        real(dp) :: elliptic

        if (any(eta < 0.0_dp)) then
            print *, "Error: eta must be non-negative!"
            error stop
        end if
        if (eps_0 >= 0.0_dp) then
            print *, "Error: eps_0 must be negative for the analytic tau calculation!"
            print *, "eps_0 = ", eps_0
            error stop
        end if

        k_squared = (1.0_dp - eta*B_min)/(2.0_dp*eta*abs(B_0*eps_0))
        if (any(k_squared < -machine_eps) .or. any(k_squared > 1.0_dp)) then
            print *, "Error: k_squared must be between 0 and 1!"
            print *, "k_squared = ", k_squared
            print *, "eta = ", eta
            print *, "eps_0 = ", eps_0
            print *, "B_0 = ", B_0
            print *, "N_tor = ", N_tor
            print *, "B_min = ", B_min
            error stop
        end if

        do idx = 1, size(eta)
            elliptic = calc_elliptic(k_squared(idx))
            ! tau = bounce_time*v/(B_varphi + iota*B_vartheta)
            tau(idx) = 8.0_dp*elliptic/sqrt(2.0_dp*abs(B_0*eps_0)*eta(idx))/(B_0*N_tor)
        end do

    end function calc_tau

    function calc_elliptic(k_squared, alpha_squared) result(elliptic)
        use constants, only: machine_eps
        real(dp), intent(in) :: k_squared
        real(dp), intent(in), optional :: alpha_squared
        real(dp) :: elliptic

        real(dp) :: agm, alpha_sq
        real(dp), allocatable :: Q_sequence(:)

        if (k_squared < -machine_eps .or. k_squared > 1.0_dp) then
            print *, "Error: k_squared must be between 0 and 1!"
            print *, "k_squared = ", k_squared
            error stop
        end if
        alpha_sq = 0.0_dp
        if (present(alpha_squared)) alpha_sq = alpha_squared
        if (alpha_sq >= 1.0_dp) then
            print *, "Error: alpha_squared must be smaller than 1!"
            print *, "alpha_squared = ", alpha_sq
            error stop
        end if

        ! Elliptic integral of the third kind via AGM (DLMF 19.8.7):
        ! integral_0^pi/2 dx/((1 - alpha^2*sin^2(x))*sqrt(1 - k^2*sin^2(x)))
        !     = pi/(4*AGM(1, sqrt(1 - k^2)))*(2 + alpha^2/(1 - alpha^2)*sum(Q_n))
        call calc_agm_and_Q_sequence(1.0_dp, sqrt(1.0_dp - k_squared), &
                                     sqrt(1.0_dp - alpha_sq), 1.0_dp, &
                                     agm, Q_sequence)
        elliptic = pi/(4.0_dp*agm) &
                   *(2.0_dp + alpha_sq/(1.0_dp - alpha_sq)*sum(Q_sequence))

    end function calc_elliptic

    !> Compute the arithmetic-geometric mean (AGM)
    !> as well as Q sidesequence (DLMF 19.8.8):
    !> eps_n = (p_n^2 - a_n*b_n)/(p_n^2 + a_n*b_n)
    !> Q_(n+1) = Q_n*eps_n/2, p_(n+1) = (p_n^2 + a_n*b_n)/(2*p_n)
    subroutine calc_agm_and_Q_sequence(a, b, p, Q, agm, Q_sequence, reltol)
        real(dp), intent(in) :: a, b, p, Q
        real(dp), intent(out) :: agm
        real(dp), dimension(:), allocatable, intent(out) :: Q_sequence
        real(dp), intent(in), optional :: reltol

        real(dp) :: a_n, b_n, a_next, b_next
        real(dp) :: p_n, p_next, eps_n
        integer :: n
        real(dp) :: tol
        integer, parameter :: max_iter = 1000
        real(dp), dimension(max_iter) :: Q_temp

        if (present(reltol)) then
            tol = reltol*(abs(a) + abs(b))
        else
            tol = 1e-12_dp
        end if

        a_n = a
        b_n = b
        p_n = p
        Q_temp(1) = Q

        n = 1

        do while (abs(a_n - b_n) > tol)
            eps_n = (p_n**2.0_dp - a_n*b_n)/(p_n**2.0_dp + a_n*b_n)
            Q_temp(n + 1) = 0.5_dp*Q_temp(n)*eps_n
            p_next = (p_n**2.0_dp + a_n*b_n)/(2.0_dp*p_n)
            a_next = 0.5_dp*(a_n + b_n)
            b_next = sqrt(a_n*b_n)
            a_n = a_next
            b_n = b_next
            p_n = p_next
            n = n + 1
            if (n > max_iter) then
                print *, "Error: AGM did not converge!"
                error stop
            end if
        end do

        agm = a_n
        allocate (Q_sequence(n))
        Q_sequence = Q_temp(1:n)

    end subroutine calc_agm_and_Q_sequence

    function calc_deeply_trapped_tau(eta, eps_0, B_0, N_tor) result(tau)
        real(dp), dimension(:), intent(in) :: eta
        real(dp), intent(in) :: eps_0, B_0, N_tor

        real(dp), dimension(size(eta)) :: tau
        real(dp), dimension(size(eta)) :: k_squared
        real(dp) :: d2B_dphi2, B_min

        d2B_dphi2 = abs(B_0*eps_0*N_tor**2.0_dp)
        B_min = B_0*(1.0_dp - abs(eps_0))
        k_squared = (1.0_dp - B_min*eta)/(2.0_dp*eta*abs(B_0*eps_0))
        tau = 4.0_dp*pi*sqrt(0.5_dp/(d2B_dphi2*eta))*(1.0_dp + k_squared/4.0_dp)/B_0

    end function calc_deeply_trapped_tau

end program test_precession_analytic
