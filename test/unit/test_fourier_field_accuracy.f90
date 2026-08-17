program plot_fourier_field_accuracy
    use myplot_module, only: myplot
    use constants, only: dp, pi
    use horner_fourier_field, only: horner_fourier_field_t
    use fourier_field, only: fourier_field_t, fourier_field_init

    implicit none

    integer, parameter :: n_m = 6, n_n = 3, n_total = n_m*n_n
    integer :: m(n_total), n(n_total)
    real(dp) :: B_modes(n_total)
    integer :: m_idx, n_idx, k_mode

    real(dp), parameter :: B_vartheta_covariant = 0.0_dp
    real(dp), parameter :: B_varphi_covariant = 1.0_dp
    integer, parameter :: n_sweep = 8
    integer :: n_grids(n_sweep)
    real(dp) :: n_grids_dp(n_sweep)

    real(dp) :: max_rel_err_B(n_sweep)
    real(dp) :: max_rel_err_dB(n_sweep)

    integer, parameter :: n_check = 100
    real(dp) :: theta_check(n_check), phi_check(n_check)
    real(dp) :: B_ref(n_check), dB_ref(3, n_check)
    real(dp) :: B_spl, dB_spl(3)
    real(dp) :: rnd, rel_err_B, rel_err_dB
    real(dp) :: ref_line_B(n_sweep), ref_line_dB(n_sweep)
    integer, parameter :: nfp = 2

    logical :: test_failed
    real(dp) :: tol = 1e-1
    integer :: n_seed
    integer, allocatable :: seed(:)

    type(horner_fourier_field_t) :: hfield
    type(fourier_field_t) :: sfield
    integer :: i, k

    logical :: should_plot = .false.
    type(myplot) :: plt

    test_failed = .false.

    n_grids = [8, 16, 32, 64, 128, 256, 512, 1024] + 1

    k_mode = 0
    do m_idx = 1, n_m
        do n_idx = 1, n_n
            k_mode = k_mode + 1
            m(k_mode) = m_idx - 1
            n(k_mode) = n_idx - 2         ! n = -1, 0, 1
            if (m(k_mode) == 0 .and. n(k_mode) == 0) then
                B_modes(k_mode) = 1.0_dp
            else
                B_modes(k_mode) = 0.01_dp
            end if
        end do
    end do

    call hfield%horner_fourier_field_init(m, n, B_modes, nfp=nfp)

    ! for reproducibile test results we set the seed
    call random_seed(size=n_seed)
    allocate (seed(n_seed))
    seed = 123456789
    call random_seed(put=seed)
    do i = 1, n_check
        call random_number(rnd); theta_check(i) = rnd*2.0_dp*pi
        call random_number(rnd); phi_check(i) = rnd*2.0_dp*pi/nfp
    end do

    do i = 1, n_check
        call hfield%compute_B_and_dB_dx(theta_check(i), phi_check(i), &
                                        B_ref(i), dB_ref(:, i))
    end do

    do k = 1, n_sweep
        n_grids_dp(k) = real(n_grids(k), dp)

        call fourier_field_init(sfield, m, n, B_modes, &
                                B_vartheta_covariant, B_varphi_covariant, &
                                nfp=nfp, n_grid_in=n_grids(k))

        max_rel_err_B(k) = 0.0_dp
        max_rel_err_dB(k) = 0.0_dp

        do i = 1, n_check
            call sfield%compute_B_and_dB_dx(theta_check(i), phi_check(i), &
                                            B_spl, dB_spl)
            rel_err_B = abs(B_spl - B_ref(i))/abs(B_ref(i))
            rel_err_dB = maxval(abs(dB_spl(2:3) - dB_ref(2:3, i)))/ &
                         maxval(abs(dB_ref(2:3, i)))
            max_rel_err_B(k) = max(max_rel_err_B(k), rel_err_B)
            max_rel_err_dB(k) = max(max_rel_err_dB(k), rel_err_dB)
        end do

        print *, "n_grid =", n_grids(k), &
            "  rel_err_B =", max_rel_err_B(k), &
            "  rel_err_dB =", max_rel_err_dB(k)
    end do

    ! --- Reference lines n^{-6}, anchored to last (well-converged) point ---
    ref_line_B = max_rel_err_B(n_sweep)*(n_grids_dp(n_sweep)/n_grids_dp)**6
    ref_line_dB = max_rel_err_dB(n_sweep)*(n_grids_dp(n_sweep)/n_grids_dp)**5

    do k = 1, n_sweep
        if (n_grids(k) >= 30) then ! need sufficient points for asymptotics
            if (max_rel_err_B(k) > (1.0_dp + tol)*ref_line_B(k) &
                .or. max_rel_err_dB(k) > (1.0_dp + tol)*ref_line_dB(k)) then
                print *, "test failed: fourier field accuracy for n_grid = ", n_grids(k)
                print *, " rel_err_B = ", max_rel_err_B(k)
                print *, " expected rel_err_B = ", ref_line_B(k)
                print *, " rel_err_dB = ", max_rel_err_dB(k)
                print *, " expected rel_err_dB = ", ref_line_dB(k)
                test_failed = .true.
            end if
        end if
    end do

    if (should_plot) then
        call plt%initialize( &
            xlabel="n grid", &
            ylabel="max relative error", &
            title="fourier_field_t accuracy vs. grid size (order-5 spline)", &
            legend=.true.)
        call plt%add_plot(n_grids_dp, max_rel_err_B, &
                          label="B", &
                          linestyle="-o", &
                          xscale="log", &
                          yscale="log")
        call plt%add_plot(n_grids_dp, max_rel_err_dB, &
                          label="dB/dx", &
                          linestyle="-s")
        call plt%add_plot(n_grids_dp, ref_line_B, &
                          label="n^-6 (B)", &
                          linestyle="--")
        call plt%add_plot(n_grids_dp, ref_line_dB, &
                          label="n^-5 (dB/dx)", &
                          linestyle=":")
        call plt%show()
    end if

    if (test_failed) error stop "test failed: fourier field accuracy"

end program plot_fourier_field_accuracy
