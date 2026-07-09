program plot_fourier_field_timing
    use myplot_module, only: myplot
    use constants, only: dp, pi
    use boozer_field, only: boozer_field_t
    use horner_fourier_field, only: horner_fourier_field_t
    use fourier_field, only: fourier_field_t, fourier_field_init

    implicit none

    character(len=*), parameter :: nc_filename = &
                      "input/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc"
    type(boozer_field_t) :: bfield
    type(horner_fourier_field_t) :: hfield
    type(fourier_field_t) :: sfield
    real(dp), parameter :: stor = 0.5_dp

    ! 50x50 fake Fourier modes: B = sum_{m,n} B_mn * cos(m*theta - n*phi)
    ! m = 0..49, n = -24..25, stored as flat lists of length n_m*n_n
    integer, parameter :: n_m = 50, n_n = 50, n_total = n_m*n_n
    integer :: m(n_total), n(n_total)
    real(dp) :: B_modes(n_total)
    integer :: m_idx, n_idx, k_mode
    real(dp), parameter :: B_vartheta_covariant = 0.0_dp
    real(dp), parameter :: B_varphi_covariant = 1.0_dp

    integer, parameter :: n_trials = 4
    integer :: eval_counts(n_trials)
    real(dp) :: eval_counts_dp(n_trials)
    real(dp) :: t_boozer(n_trials), t_horner(n_trials), t_spline(n_trials)

    real(dp) :: theta, phi, B_mod, dB_dx(3)
    real(dp) :: t_start, t_end
    integer :: k, i, eval
    real(dp) :: rnd
    integer :: j, tmp_m, tmp_n
    real(dp) :: tmp_B

    type(myplot) :: plt

    eval_counts = [100, 1000, 10000, 100000]

    ! Dominant B_00 = 1, others decay as 1/(m+|n|+2)
    k_mode = 0
    do m_idx = 1, n_m
        do n_idx = 1, n_n
            k_mode = k_mode + 1
            m(k_mode) = m_idx - 1
            n(k_mode) = n_idx - 1 - n_n/2
            if (m(k_mode) == 0 .and. n(k_mode) == 0) then
                B_modes(k_mode) = 1.0_dp
            else
                B_modes(k_mode) = 0.01_dp/ &
                                  real(m(k_mode) + abs(n(k_mode)) + 2, dp)
            end if
        end do
    end do

    ! --- Shuffle mode lists (Fisher-Yates) to test sorting in fourier_field_init ---
    do k_mode = n_total, 2, -1
        call random_number(rnd)
        j = int(rnd*k_mode) + 1
        tmp_m = m(k_mode); m(k_mode) = m(j); m(j) = tmp_m
        tmp_n = n(k_mode); n(k_mode) = n(j); n(j) = tmp_n
        tmp_B = B_modes(k_mode); B_modes(k_mode) = B_modes(j); B_modes(j) = tmp_B
    end do

    call bfield%init_from_vmec(nc_filename, &
                                  radial_spline_order=5, &
                                  angular_spline_order=5, &
                                  grid_refinement=3)
    call bfield%fix_to_surface(stor)

    call hfield%horner_fourier_field_init(m, n, B_modes)

    call fourier_field_init(sfield, m, n, B_modes, &
                            B_vartheta_covariant, B_varphi_covariant)

    do k = 1, n_trials
        eval = eval_counts(k)
        eval_counts_dp(k) = real(eval, dp)

        ! Time boozer_field_t
        call cpu_time(t_start)
        do i = 1, eval
            theta = real(i, dp)/real(eval, dp)*2.0_dp*pi
            phi = real(i, dp)/real(eval, dp)*2.0_dp*pi
            call bfield%compute_B_and_dB_dx(theta, phi, B_mod, dB_dx)
        end do
        call cpu_time(t_end)
        t_boozer(k) = t_end - t_start

        ! Time horner_fourier_field_t
        call cpu_time(t_start)
        do i = 1, eval
            theta = real(i, dp)/real(eval, dp)*2.0_dp*pi
            phi = real(i, dp)/real(eval, dp)*2.0_dp*pi
            call hfield%compute_B_and_dB_dx(theta, phi, B_mod, dB_dx)
        end do
        call cpu_time(t_end)
        t_horner(k) = t_end - t_start

        ! Time fourier_field_t (spline)
        call cpu_time(t_start)
        do i = 1, eval
            theta = real(i, dp)/real(eval, dp)*2.0_dp*pi
            phi = real(i, dp)/real(eval, dp)*2.0_dp*pi
            call sfield%compute_B_and_dB_dx(theta, phi, B_mod, dB_dx)
        end do
        call cpu_time(t_end)
        t_spline(k) = t_end - t_start
    end do

    call plt%initialize( &
        xlabel="Number of evaluations", &
        ylabel="CPU time (s)", &
        title="compute_B_and_dB_dx timing", &
        legend=.true.)
    call plt%add_plot(eval_counts_dp, t_boozer, &
                      label="boozer_field_t (QH)", &
                      linestyle="-o", &
                      xscale="log", &
                      yscale="log")
    call plt%add_plot(eval_counts_dp, t_horner, &
                      label="horner_fourier_field_t (50x50)", &
                      linestyle="-s")
    call plt%add_plot(eval_counts_dp, t_spline, &
                      label="fourier_field_t spline (50x50)", &
                      linestyle="-^")
    call plt%show()

end program plot_fourier_field_timing
