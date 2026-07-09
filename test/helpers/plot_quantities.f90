module plot_quantities
    use constants, only: dp, pi
    use myplot_module, only: myplot
    use fieldline_mod, only: fieldline_t, flock_of_fieldlines_t
    use utils, only: linspace
    use field_base, only: field_t
    use fieldline_labels, only: fourier_transform_over_label
    use fieldline_labels, only: fieldline_modes_t
    use fieldline_labels, only: modes_t

    implicit none

    type external_data_t
        real(dp), dimension(:), allocatable :: x
        real(dp), dimension(:), allocatable :: y
        character(len=1024) :: label
    end type external_data_t

contains

    subroutine plot_maxima_over_label(fieldlines, iota_p)
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp), intent(in) :: iota_p

        type(myplot) :: plt
        real(dp), dimension(size(fieldlines)) :: xi_0_shifted

        call plt%initialize(xlabel="$\xi_0 - \iota_p$ [$\pi$]", &
                            ylabel="$B$ [T]", &
                            figsize=[7, 7], &
                            legend=.true.)
        xi_0_shifted = modulo(fieldlines%xi_0 - iota_p, 2.0_dp*pi)
        call plt%add_plot(xi_0_shifted/pi, &
                          fieldlines%B_max(1), &
                          label="left $B_{max}$", &
                          linestyle="-o")
        call plt%add_plot(xi_0_shifted/pi, &
                          fieldlines%B_max(2), &
                          label="right $B_{max}$", &
                          linestyle="--o")
        call plt%show()
    end subroutine plot_maxima_over_label

    subroutine plot_field_along_chi_line(field, chi, M_pol, N_tor)
        use utils, only: linspace
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: chi, M_pol, N_tor

        type(myplot) :: plt
        integer, parameter :: n_points = 100
        real(dp), dimension(:), allocatable :: theta, phi, B

        character(len=1024) :: label
        integer :: current

        allocate (theta(n_points), phi(n_points), B(n_points))

        call linspace(0.0_dp, 2.0_dp*pi, n_points, theta)
        phi = (theta*M_pol - chi)/N_tor

        do current = 1, n_points
            call field%compute_B_mod(theta(current), phi(current), B(current))
        end do

        call plt%initialize(xlabel="$\vartheta$ [$\pi$]", &
                            ylabel="$B$ [T]", &
                            figsize=[7, 7], &
                            legend=.true.)
        write (label, "(A17,ES10.3E2)") "$B$ along $\chi=$", chi
        call plt%add_plot(theta/pi, B, label=label, linestyle="b-")
        call plt%show()

        deallocate (theta, phi, B)
    end subroutine plot_field_along_chi_line

    subroutine plot_spectra(flock)
        type(flock_of_fieldlines_t), intent(in) :: flock

        type(fieldline_modes_t) :: modes

        integer :: n_modes
        real(dp), dimension(:), allocatable :: deviation_B
        integer :: current
        real(dp) :: iota_p

        type(myplot) :: plt
        character(len=1024) :: label

        call fourier_transform_over_label(flock, modes)
        n_modes = size(modes%delta_eta%cos_coeffs)

        call plt%initialize(xlabel="mode", &
                            ylabel="amplitude", &
                            legend=.true.)

        call plt%add_plot(modes%delta_aspect_ratio%mode_numbers, &
                          modes%delta_aspect_ratio%cos_coeffs, &
                          label="cos $\Delta A$", &
                          linestyle="-ro")
        call plt%add_plot(modes%delta_eta%mode_numbers, &
                          modes%delta_eta%cos_coeffs, &
                          label="cos $\Delta \eta$", &
                          linestyle="-bo")
        call plt%add_plot(modes%radial_drift%mode_numbers, &
                          modes%radial_drift%sin_coeffs, &
                          label="sin $I_\mathcal{V}$", &
                          linestyle="-ko")

        call plt%show()

    end subroutine plot_spectra

    subroutine plot_deviation_spectrum(flock)
        use fit_functions, only: S_B
        type(flock_of_fieldlines_t), intent(in) :: flock

        type(fieldline_modes_t) :: modes

        integer :: n_modes
        real(dp), dimension(:), allocatable :: deviation_B
        integer :: current
        real(dp) :: iota_p

        type(myplot) :: plt
        character(len=1024) :: label

        call fourier_transform_over_label(flock, modes)
        iota_p = flock%iota_p
        n_modes = size(modes%delta_eta%cos_coeffs)

        allocate (deviation_B(n_modes))

        call plt%initialize(xlabel="mode", &
                            ylabel="amplitude", &
                            legend=.true.)

        do current = 1, n_modes
            deviation_B(current) = modes%radial_drift%sin_coeffs(current)* &
                                   modes%delta_eta%cos_coeffs(current)* &
                                   S_B(iota_p*modes%delta_eta%mode_numbers(current))
        end do

        call plt%add_plot(modes%delta_eta%mode_numbers, &
                          deviation_B, &
                          label="spectrum of deviation due to $\Delta \eta$", &
                          linestyle="-o")

        call plt%show()

        deallocate (deviation_B)
    end subroutine plot_deviation_spectrum

    subroutine plot_delta_eta_modes(flock)
        type(flock_of_fieldlines_t), intent(in) :: flock

        type(fieldline_modes_t) :: fieldline_modes

        integer, parameter :: n_points = 100
        real(dp), dimension(:), allocatable :: xi_mid, xi, delta_eta
        integer :: max_mode
        integer, parameter :: skip = 3

        type(myplot) :: plt
        character(len=1024) :: label

        call fourier_transform_over_label(flock, fieldline_modes)

        allocate (xi_mid(n_points), xi(n_points), delta_eta(n_points))
        call linspace(0.0_dp, 2.0_dp*pi, n_points, xi_mid)
        xi = xi_mid - flock%iota_p

        call plt%initialize(xlabel="$\xi_{mid}[\pi]$", &
                            ylabel="$\Delta \eta$", &
                            legend=.true.)

        do max_mode = 1, size(fieldline_modes%delta_eta%cos_coeffs), skip
            delta_eta = eval_modes(xi, fieldline_modes%delta_eta, max_mode)
            write (label, "(I3,A6)") max_mode, " modes"
            call plt%add_plot(xi_mid/pi, &
                              delta_eta, &
                              label=label, &
                              linestyle="-")
        end do

        call plt%add_plot(flock%fieldlines%xi_0/pi, &
                          flock%fieldlines%delta_eta, &
                          label="original", &
                          linestyle="ko")

        call plt%show()

        deallocate (xi_mid, delta_eta)
    end subroutine plot_delta_eta_modes

    elemental function eval_modes(x, modes, max_mode)
        real(dp), intent(in) :: x
        type(modes_t), intent(in) :: modes
        integer, intent(in) :: max_mode

        real(dp) :: eval_modes

        eval_modes = sum(modes%cos_coeffs(1:max_mode) &
                         *cos(modes%mode_numbers(1:max_mode)*x)) + &
                     sum(modes%sin_coeffs(1:max_mode) &
                         *sin(modes%mode_numbers(1:max_mode)*x))
    end function eval_modes

    subroutine plot_fieldlines_over_field(fieldlines, field)

        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        class(field_t), intent(in) :: field

        type(myplot) :: plt
        integer :: current
        integer :: n_fieldlines
        real(dp), dimension(2) :: phi_max, theta_max
        real(dp), dimension(4) :: phi_limits, theta_limits
        integer, parameter :: n_points = 300
        real(dp), dimension(n_points) :: phi, theta
        real(dp), dimension(:, :), allocatable :: B_mesh
        integer :: theta_idx, phi_idx
        character(len=100) :: label

        allocate (B_mesh(n_points, n_points))

        n_fieldlines = size(fieldlines)

        call plt%initialize(xlabel="$\varphi$", ylabel="$\vartheta$")

        phi_limits(1) = fieldlines(1)%phi_max(1)
        phi_limits(2) = fieldlines(1)%phi_max(2)
        phi_limits(3) = fieldlines(n_fieldlines)%phi_max(1)
        phi_limits(4) = fieldlines(n_fieldlines)%phi_max(2)

        theta_limits(1:2) = fieldlines(1)%get_theta(fieldlines(1)%phi_max)
        theta_limits(3:4) = fieldlines(n_fieldlines)%get_theta( &
                            fieldlines(n_fieldlines)%phi_max)

        do current = 1, size(fieldlines)
            phi_max = fieldlines(current)%phi_max
            call linspace(phi_max(1), &
                          phi_max(2), &
                          n_points, &
                          phi)
            theta = fieldlines(current)%get_theta(phi)
            write (label, '(F0.1)') fieldlines(current)%xi_0
            call plt%add_plot(phi, theta, &
                              label=label, &
                              linestyle="k.", &
                              linewidth=1)
            theta_max = fieldlines(current)%get_theta(phi_max)
            call plt%add_plot(phi_max, theta_max, &
                              label=label, &
                              linestyle="ro", &
                              markersize=5)
        end do

        call linspace(minval(phi_limits), maxval(phi_limits), n_points, phi)
        call linspace(minval(theta_limits), maxval(theta_limits), n_points, theta)
        do theta_idx = 1, n_points
            do phi_idx = 1, n_points
                call field%compute_B_mod(theta(theta_idx), &
                                         phi(phi_idx), &
                                         B_mesh(theta_idx, phi_idx))
            end do
        end do
        call plt%add_contour(phi, theta, transpose(B_mesh), &
                             levels=20, &
                             colorbar=.true., &
                             filled=.true.)
        call plt%show()

    end subroutine plot_fieldlines_over_field

    subroutine plot_fieldlines_over_field_chi_xi(fieldlines, field, M_pol, N_tor, &
                                                 nfp)

        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: M_pol, N_tor, nfp

        type(myplot) :: plt
        integer :: current
        integer :: n_fieldlines
        real(dp), dimension(2) :: phi_max, theta_max
        real(dp), dimension(2) :: chi_max, xi_max
        real(dp), dimension(4) :: phi_limits, theta_limits
        real(dp), dimension(4) :: chi_limits, xi_limits
        integer, parameter :: n_points = 300, linewidth = 2
        real(dp), dimension(n_points) :: phi, theta
        real(dp), dimension(n_points) :: chi, xi
        real(dp) :: phi_temp, theta_temp
        real(dp), dimension(:, :), allocatable :: B_mesh
        integer :: xi_idx, chi_idx
        character(len=100) :: label

        allocate (B_mesh(n_points, n_points))

        n_fieldlines = size(fieldlines)

        call plt%initialize(xlabel="$\chi$", ylabel="$\xi$")

        phi_limits(1) = fieldlines(1)%phi_max(1)
        phi_limits(2) = fieldlines(1)%phi_max(2)
        phi_limits(3) = fieldlines(n_fieldlines)%phi_max(1)
        phi_limits(4) = fieldlines(n_fieldlines)%phi_max(2)

        theta_limits(1:2) = fieldlines(1)%get_theta(fieldlines(1)%phi_max)
        theta_limits(3:4) = fieldlines(n_fieldlines)%get_theta( &
                            fieldlines(n_fieldlines)%phi_max)

        call convert_to_chi_xi(theta_limits, phi_limits, &
                               M_pol, N_tor, nfp, &
                               xi_limits, chi_limits)
        do current = 1, size(fieldlines)
            phi_max = fieldlines(current)%phi_max
            call linspace(phi_max(1), &
                          phi_max(2), &
                          n_points, &
                          phi)
            theta = fieldlines(current)%get_theta(phi)
            call convert_to_chi_xi(theta, phi, M_pol, N_tor, nfp, xi, chi)
            write (label, '(F0.1)') fieldlines(current)%xi_0
            call plt%add_plot(chi, xi, &
                              label=label, &
                              linestyle="k.", &
                              linewidth=1)
            theta_max = fieldlines(current)%get_theta(phi_max)
            call convert_to_chi_xi(theta_max, phi_max, M_pol, N_tor, nfp, xi_max, &
                                   chi_max)
            call plt%add_plot(chi_max, xi_max, &
                              label=label, &
                              linestyle="ro", &
                              markersize=5)
        end do

        call linspace(minval(chi_limits), maxval(chi_limits), n_points, chi)
        call linspace(minval(xi_limits), maxval(xi_limits), n_points, xi)
        do xi_idx = 1, n_points
            do chi_idx = 1, n_points
                call convert_to_theta_phi(xi(xi_idx), chi(chi_idx), M_pol, N_tor, nfp, &
                                          theta_temp, phi_temp)
                call field%compute_B_mod(theta_temp, &
                                         phi_temp, &
                                         B_mesh(xi_idx, chi_idx))
            end do
        end do
        call plt%add_contour(chi, xi, transpose(B_mesh), &
                             levels=20, &
                             colorbar=.true., &
                             filled=.true.)

        call plt%add_plot([chi(1), chi(n_points)], [1.0_dp, 1.0_dp]*M_pol/nfp*pi, &
                          label="", &
                          linestyle="k-", &
                          linewidth=linewidth)
        call plt%add_plot(-pi*[1.0_dp, 1.0_dp], [xi(1), xi(n_points)], &
                          label="", &
                          linestyle="k-", &
                          linewidth=linewidth)
        call plt%add_plot([chi(1), chi(n_points)], -[1.0_dp, 1.0_dp]*M_pol/nfp*pi, &
                          label="", &
                          linestyle="k-", &
                          linewidth=linewidth)
        call plt%add_plot(pi*[1.0_dp, 1.0_dp], [xi(1), xi(n_points)], &
                          label="", &
                          linestyle="k-", &
                          linewidth=linewidth)
        call plt%show()

    end subroutine plot_fieldlines_over_field_chi_xi

    subroutine convert_to_chi_xi(theta, phi, M_pol, N_tor, nfp, xi, chi)
        real(dp), dimension(:), intent(in) :: theta, phi
        real(dp), intent(in) :: M_pol, N_tor, nfp
        real(dp), dimension(:), intent(out) :: xi, chi

        real(dp) :: normalization

        normalization = M_pol**2.0_dp + N_tor**2.0_dp

        xi = (N_tor*theta + M_pol*phi)/normalization*nfp
        chi = M_pol*theta - N_tor*phi
    end subroutine convert_to_chi_xi

    subroutine convert_to_theta_phi(xi, chi, M_pol, N_tor, nfp, theta, phi)
        real(dp), intent(in) :: xi, chi
        real(dp), intent(in) :: M_pol, N_tor, nfp
        real(dp), intent(out) :: theta, phi

        real(dp) :: normalization

        normalization = M_pol**2.0_dp + N_tor**2.0_dp

        theta = N_tor*xi/nfp + M_pol*chi/normalization
        phi = M_pol*xi/nfp - N_tor*chi/normalization
    end subroutine convert_to_theta_phi

    subroutine plot_phi_max_over_xi_0(fieldlines, M_pol, nfp)
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp), intent(in) :: M_pol, nfp
        type(myplot) :: plt

        real(dp), dimension(size(fieldlines)) :: xi_0, phi_l, phi_r

        xi_0 = fieldlines%xi_0
        phi_l = fieldlines%phi_max(1)
        phi_r = fieldlines%phi_max(2)

        call plt%initialize(xlabel="$\xi_0$", &
                            ylabel="$\varphi_\mathrm{max}$", &
                            legend=.true.)

        call plt%add_plot(xi_0, phi_l, &
                          label="$\varphi_l$", &
                          linestyle="bo-", &
                          linewidth=1)
        call plt%add_plot(xi_0, phi_r, &
                          label="$\varphi_r$", &
                          linestyle="ro-", &
                          linewidth=1)
        call plt%add_plot(xi_0, phi_l - M_pol/nfp*xi_0, &
                          label="$\varphi_l - M/N_p \xi_0$", &
                          linestyle="b--", &
                          linewidth=1)
        call plt%add_plot(xi_0, phi_r - M_pol/nfp*xi_0, &
                          label="$\varphi_r - M/N_p \xi_0$", &
                          linestyle="r--", &
                          linewidth=1)

        call plt%show()

    end subroutine plot_phi_max_over_xi_0

    subroutine plot_fieldline_over_local_drift(fieldline, field, eta, interval)
        use fieldline_integrands, only: local_radial_drift
        type(fieldline_t), intent(in) :: fieldline
        class(field_t), intent(in) :: field
        real(dp), intent(in) :: eta
        real(dp), dimension(2), intent(in) :: interval

        type(myplot) :: plt
        integer, parameter :: n_points = 300
        real(dp), dimension(2) :: phi_limits, theta_limits
        real(dp) :: phi_range
        real(dp), dimension(2) :: phi_ends, theta_ends
        real(dp), dimension(n_points) :: phi, theta
        real(dp), dimension(:, :), allocatable :: drift_mesh
        integer :: theta_idx, phi_idx
        character(len=100) :: label, title
        character(len=*), parameter :: cmap = "coolwarm"

        allocate (drift_mesh(n_points, n_points))

        write (title, "(A,F4.2)") "Local radial drift for $\eta=$", eta

        call plt%initialize(xlabel="$\varphi$", ylabel="$\vartheta$", &
                            title=title, fontsize=20)

        phi_limits = interval
        theta_limits(1) = -pi
        theta_limits(2) = pi
        phi_range = abs(phi_limits(2) - phi_limits(1))

        phi_ends(1) = fieldline%phi_0 - 0.4_dp*phi_range
        phi_ends(2) = fieldline%phi_0 + 0.4_dp*phi_range

        call linspace(phi_ends(1), phi_ends(2), n_points, phi)
        theta = fieldline%get_theta(phi)
        write (label, '(F0.1)') fieldline%xi_0
        call plt%add_plot(phi, theta, &
                          label=label, &
                          linestyle="k.", &
                          linewidth=1)
        theta_ends = fieldline%get_theta(phi_ends)
        call plt%add_plot(phi_ends, theta_ends, &
                          label=label, &
                          linestyle="ro", &
                          markersize=10)

        call linspace(minval(phi_limits), maxval(phi_limits), n_points, phi)
        call linspace(minval(theta_limits), maxval(theta_limits), n_points, theta)
        do theta_idx = 1, n_points
            do phi_idx = 1, n_points
                drift_mesh(theta_idx, phi_idx) = local_radial_drift(field, &
                                                                    theta(theta_idx), &
                                                                    phi(phi_idx), &
                                                                    eta)
            end do
        end do
        call plt%add_contour(phi, theta, transpose(drift_mesh), &
                             levels=20, &
                             colorbar=.true., &
                             filled=.true., &
                             cmap=cmap)
        call plt%show()

    end subroutine plot_fieldline_over_local_drift

    subroutine plot_local_drift_over_fieldline(field, &
                                               fieldline, &
                                               etas, &
                                               interval)
        use fieldline_integrands, only: local_radial_drift
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        class(field_t), intent(in) :: field
        type(fieldline_t), intent(in) :: fieldline
        real(dp), dimension(:), intent(in) :: etas
        real(dp), intent(in), optional :: interval(2)

        integer, parameter :: n_points = 1001
        character(len=*), parameter :: cmap = "coolwarm"
        integer :: n_etas
        real(dp), dimension(:), allocatable :: phi, theta, B, eta_level
        real(dp), dimension(:, :), allocatable :: drift
        character(len=1024) :: label
        character(len=1024) :: title
        integer :: current, id_eta

        type(myplot) :: plt

        allocate (phi(n_points), theta(n_points), B(n_points))
        n_etas = size(etas)
        allocate (eta_level(n_points))
        allocate (drift(n_points, n_etas))

        if (present(interval)) then
            call linspace(interval(1), interval(2), n_points, phi)
        else
            call linspace(0.0_dp, 2.0_dp*pi, n_points, phi)
        end if

        theta = fieldline%get_theta(phi)

        write (title, "(A,A,F5.2,A)") "Radial drift velocity ", &
            " along fieldline $\xi_0=$", fieldline%xi_0/pi, "[$\pi$]"
        call plt%initialize(xlabel="$\varphi$", &
                            ylabel="$1/\eta$ [T]", &
                            legend=.true., &
                            figsize=[10, 10], &
                            title=title, &
                            fontsize=20)

        do current = 1, n_points
            call field%compute_B_mod(theta(current), phi(current), B(current))
            do id_eta = 1, n_etas
                if (B(current)*etas(id_eta) > 1.0_dp) then
                    ! set to nan with intrinsic
                    drift(current, id_eta) = ieee_value(1.0_dp, ieee_quiet_nan)
                else
                    drift(current, id_eta) = local_radial_drift(field, &
                                                                theta(current), &
                                                                phi(current), &
                                                                etas(id_eta))
                end if
            end do
        end do

        eta_level = etas(1)
        call plt%add_colored_line(phi, 1.0_dp/eta_level, drift(:, 1), &
                                  cmap=cmap, &
                                  clabel="$\mathcal{V}$", linewidth=3)
        do id_eta = 2, n_etas
            eta_level = etas(id_eta)
            call plt%add_colored_line(phi, 1.0_dp/eta_level, drift(:, id_eta), &
                                      cmap=cmap, &
                                      linewidth=3)
        end do

        write (label, "(A,F6.2,A)") "$B$ [T]"
        call plt%add_plot(phi, B, label=label, linestyle="k-", linewidth=3)

        call plt%show()

        deallocate (phi, theta, drift, B, eta_level)

    end subroutine plot_local_drift_over_fieldline

    subroutine plot_delta_eta(fieldlines, iota_p, delta_eta_1)
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp), intent(in) :: iota_p
        real(dp), intent(in), optional :: delta_eta_1

        type(myplot) :: plt

        call plt%initialize(xlabel="$\xi_{mid}$ [$\pi$]", &
                            ylabel="$\Delta \eta$", &
                            legend=.true.)

        call plt%add_plot(fieldlines%xi_0/pi, &
                          fieldlines%delta_eta, &
                          label="$\Delta \eta$", &
                          linestyle="-")

        if (present(delta_eta_1)) then
            call plt%add_plot(fieldlines%xi_0/pi, &
                              abs(delta_eta_1) + &
                              delta_eta_1*cos(fieldlines%xi_0 - iota_p), &
                              label="$\Delta \eta$ approx analytic", &
                              linestyle="--")
        end if

        call plt%show()
    end subroutine plot_delta_eta

    subroutine plot_delta_A(fieldlines, delta_A_1_analytic)
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp), intent(in), optional :: delta_A_1_analytic

        type(myplot) :: plt

        call plt%initialize(xlabel="$\xi_{mid} [\pi]$", &
                            ylabel="$\Delta A$", &
                            legend=.true.)

        call plt%add_plot(fieldlines%xi_0/pi, &
                          fieldlines%delta_aspect_ratio, &
                          label="$\Delta A$", &
                          linestyle="-")
        if (present(delta_A_1_analytic)) then
            call plt%add_plot(fieldlines%xi_0/pi, &
                              delta_A_1_analytic*cos(fieldlines%xi_0), &
                              label="$\Delta A$ approx analytic", &
                              linestyle="-")
        end if
        call plt%show()
    end subroutine plot_delta_A

    subroutine plot_I(fieldlines, I_0_analytic, I_1_analytic, eps_0, eps_1)
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp), intent(in) :: I_0_analytic, I_1_analytic, eps_0, eps_1

        type(myplot) :: plt
        integer :: n_fieldlines
        real(dp) :: eps_ratio, I_factor, I_mean
        real(dp), dimension(size(fieldlines)) :: theta_0
        character(len=100) :: label

        n_fieldlines = size(fieldlines)
        theta_0 = fieldlines%theta_0
        eps_ratio = eps_1/abs(eps_0)
        I_mean = sum(fieldlines%integral_lambda_b_over_B_squared)/n_fieldlines
        I_factor = sum(fieldlines%integral_lambda_b_over_B_squared/ &
                       sqrt(1.0_dp - eps_ratio*cos(theta_0)))/n_fieldlines

        call plt%initialize(xlabel="$\vartheta_{mid}$", &
                            ylabel="$I$ [T$^{-2}$]", &
                            legend=.true.)
        call plt%add_plot(theta_0, &
                          fieldlines%integral_lambda_b_over_B_squared, &
                          label="$I_{numeric}$", &
                          linestyle="b-")
        write (label, "(A23,F10.8)") "$I^{mean}_{numeric} = $", I_mean
        call plt%add_plot(theta_0, &
                          I_mean*cos(0.0_dp*theta_0), &
                          label=label, &
                          linestyle="b--")
        call plt%add_plot(theta_0, &
                          I_0_analytic + I_1_analytic*cos(theta_0), &
                          label="approx $I_{analytic}$", &
                          linestyle="r-")
        call plt%add_plot(theta_0, &
                          I_0_analytic*cos(0.0_dp*theta_0), &
                          label="approx $I^{mean}_{analytic}$", &
                          linestyle="r--")
        call plt%add_plot(theta_0, &
                          I_0_analytic*sqrt(1.0_dp - eps_ratio*cos(theta_0)), &
                          label="$I_{analytic}$", &
                          linestyle="g--")
        call plt%show()

        call plt%initialize(xlabel="$\vartheta_{mid}$", &
                            ylabel="$I_{normalized}$ [1]", &
                            legend=.true.)
        call plt%add_plot(theta_0, &
                          fieldlines%integral_lambda_b_over_B_squared/I_factor &
                          - 1.0_dp, &
                          label="$I_{numeric}$", &
                          linestyle="b-")
        call plt%add_plot(theta_0, &
                          sqrt(1.0_dp - eps_ratio*cos(theta_0)) - 1.0_dp, &
                          label="$I_{analytic}$", &
                          linestyle="g--")
        call plt%show()
    end subroutine plot_I

    subroutine plot_I_factor(fieldlines, eps_0, eps_1)
        type(fieldline_t), dimension(:), intent(in) :: fieldlines
        real(dp), intent(in) :: eps_0, eps_1

        type(myplot) :: plt
        real(dp), dimension(size(fieldlines)) :: theta_0, I_factor
        integer :: n_fieldlines
        real(dp) :: I_factor_mean
        character(len=1024) :: label

        theta_0 = fieldlines%theta_0
        I_factor = fieldlines%integral_lambda_b_over_B_squared &
                   /sqrt(1.0_dp - eps_1/abs(eps_0)*cos(theta_0))
        n_fieldlines = size(fieldlines)
        I_factor_mean = sum(I_factor)/n_fieldlines

        call plt%initialize(xlabel="$\vartheta_{mid}$", &
                            ylabel="$I/\sqrt{1+\epsilon\cos{\vartheta_0}}$ [1/T$^2$]", &
                            legend=.true.)
        write (label, "(A38,F10.8)") "deviation from numeric $I_{factor} =$ ", &
            I_factor_mean
        call plt%add_plot(theta_0, &
                          I_factor/I_factor_mean, &
                          label=label, &
                          linestyle="k-")
        label = "analytic approximation $+ \mathcal{O}(\epsilon_0^2)$"
        call plt%add_plot(theta_0, &
                          -8.0_dp/(3.0_dp - 2.0_dp*eps_0)*eps_1*cos(theta_0) &
                          + 1.0_dp, &
                          label=label, &
                          linestyle="g-")
        call plt%show()
    end subroutine plot_I_factor

    subroutine plot_deviation(off_factor_A, off_factor_B, &
                              off_factor_A_analytic, off_factor_B_analytic, &
                              lambda_off_external)
        real(dp), intent(in) :: off_factor_A, off_factor_B

        real(dp), intent(in), optional :: off_factor_A_analytic
        real(dp), intent(in), optional :: off_factor_B_analytic

        type(external_data_t), intent(in), optional :: lambda_off_external

        type(myplot) :: plt
        integer, parameter :: n_points = 20
        real(dp), dimension(:), allocatable :: nu_star
        real(dp), dimension(2) :: nu_star_lim

        character(len=1024) :: label

        call plt%initialize(xlabel="$\nu_*$", &
                            ylabel="$|\lambda_{bB}|$", &
                            legend=.true., &
                            figsize=[15, 10])

        if (present(lambda_off_external)) then
            call plt%add_plot(lambda_off_external%x, &
                              lambda_off_external%y, &
                              label=lambda_off_external%label, &
                              linestyle="co", &
                              markersize=8, &
                              xscale="log", &
                              yscale="log")
            allocate (nu_star(size(lambda_off_external%x)), &
                      source=lambda_off_external%x)
        else
            allocate (nu_star(n_points))
            call linspace(0.0_dp, 8.0_dp, n_points, nu_star)
            nu_star = 0.1_dp**nu_star
        end if

        write (label, "(A35,ES10.3E2,A1)") "offset due to $\Delta A$ (factor$=$", &
            off_factor_A, ")"
        call plt%add_plot(nu_star, &
                          abs(off_factor_A)/sqrt(nu_star), &
                          label=label, &
                          linestyle="r-", &
                          xscale="log", &
                          yscale="log")

        if (present(off_factor_A_analytic)) then
            write (label, "(A31,ES10.3E2,A1)") &
                "analytic estimate (rel. diff. =", &
                abs(off_factor_A_analytic/off_factor_A - 1.0_dp), ")"
            call plt%add_plot(nu_star, &
                              abs(off_factor_A_analytic)/sqrt(nu_star), &
                              label=label, &
                              linestyle="ro", &
                              markersize=8, &
                              xscale="log", &
                              yscale="log")
        end if

        write (label, "(A38,ES10.3E2,A1)") "offset due to $\Delta \eta$ (factor$=$", &
            off_factor_B, ")"
        call plt%add_plot(nu_star, &
                          abs(off_factor_B)/nu_star, &
                          label=label, &
                          linestyle="b-", &
                          xscale="log", &
                          yscale="log")

        if (present(off_factor_B_analytic)) then
            write (label, "(A31,ES10.3E2,A1)") &
                "analytic estimate (rel. diff. =", &
                abs(off_factor_B_analytic/off_factor_B - 1.0_dp), ")"
            call plt%add_plot(nu_star, &
                              abs(off_factor_B_analytic)/nu_star, &
                              label=label, &
                              linestyle="bo", &
                              markersize=8, &
                              xscale="log", &
                              yscale="log")
        end if

        nu_star_lim(1) = minval(nu_star)
        nu_star_lim(2) = maxval(nu_star)
        call plt%add_plot(nu_star, &
                          abs(off_factor_A/sqrt(nu_star) + off_factor_B/nu_star), &
                          label="total", &
                          linestyle="c--", &
                          xscale="log", &
                          yscale="log", &
                          xlim=nu_star_lim)

        call plt%show()
    end subroutine plot_deviation

    subroutine plot_asymptotic_model(off_factor_A, off_factor_B, &
                                     lambda_SC, Lambda_S)
        real(dp), intent(in) :: off_factor_A, off_factor_B
        real(dp), intent(in) :: lambda_SC, Lambda_S

        type(myplot) :: plt
        integer, parameter :: n_points = 100
        real(dp), dimension(:), allocatable :: nu_star

        character(len=1024) :: label
        call plt%initialize(xlabel="$\nu_*$", &
                            ylabel="$|\lambda_{bB}|$", &
                            legend=.true., &
                            figsize=[15, 10])
        allocate (nu_star(n_points))
        call linspace(2.0_dp, 8.0_dp, n_points, nu_star)
        nu_star = 0.1_dp**nu_star
        write (label, "(A,ES10.3E2)") "$\Lambda_\mathrm{lm}=$", off_factor_B
        call plt%add_plot(nu_star, &
                          off_factor_B/nu_star, &
                          label=label, &
                          linestyle="None", &
                          xscale="log", &
                          yscale="linear")
        write (label, "(A,ES10.3E2)") "$\Lambda_\mathrm{bl}=$", off_factor_A
        call plt%add_plot(nu_star, &
                          off_factor_A/sqrt(nu_star), &
                          label=label, &
                          linestyle="None", &
                          xscale="log", &
                          yscale="linear")
        write (label, "(A,ES10.3E2)") "$\lambda_\mathrm{SC}=$", lambda_SC
        call plt%add_plot(nu_star, &
                          lambda_SC*nu_star/nu_star, &
                          label=label, &
                          linestyle="None", &
                          xscale="log", &
                          yscale="linear")
        write (label, "(A,ES10.3E2)") "$\Lambda_\mathrm{finite}$=", Lambda_S
        call plt%add_plot(nu_star, &
                          Lambda_S*sqrt(nu_star), &
                          label=label, &
                          linestyle="None", &
                          xscale="log", &
                          yscale="linear")
        write (label, "(A,A,A,A)") "$\lambda_\mathrm{SC} + ", &
            "\Lambda_\mathrm{finite} \sqrt{\nu_*} + ", &
            "\Lambda_\mathrm{bl}/\sqrt{\nu_*} + ", &
            "\Lambda_\mathrm{lm}/\nu_*$"
        call plt%add_plot(nu_star, &
                          lambda_SC + &
                          Lambda_S*sqrt(nu_star) + &
                          off_factor_A/sqrt(nu_star) + &
                          off_factor_B/nu_star, &
                          label=label, &
                          linestyle="k-", &
                          xscale="log", &
                          yscale="linear")
        call plt%show()
    end subroutine plot_asymptotic_model

    subroutine plot_B_along_fieldline(field, &
                                      fieldline, &
                                      interval)
        class(field_t), intent(in) :: field
        type(fieldline_t), intent(in) :: fieldline
        real(dp), intent(in), optional :: interval(2)

        integer, parameter :: n_points = 1001
        real(dp), dimension(:), allocatable :: phi, B
        character(len=1024) :: label

        type(myplot) :: plt

        allocate (phi(n_points), B(n_points))

        if (present(interval)) then
            call linspace(interval(1), interval(2), n_points, phi)
        else
            call linspace(0.0_dp, 2.0_dp*pi, n_points, phi)
        end if

        call B_mod_along_fieldline(phi, B)

        call plt%initialize(xlabel="$\varphi [\pi]$", &
                            ylabel="$B$ [T]", &
                            legend=.true., &
                            figsize=[10, 10])

        write (label, "(A25,ES10.3E2,A7)") "$\xi_\mathrm{mid}=$", &
            fieldline%xi_0/pi, "[$\pi$]"

        call plt%add_plot(phi/pi, &
                          B, &
                          label=label, &
                          linestyle="ko-", &
                          markersize=1)

        call plt%show()

        deallocate (phi, B)

    contains
        subroutine B_mod_along_fieldline(phi, B_mod)
            real(dp), dimension(:), intent(in) :: phi
            real(dp), dimension(:), intent(out) :: B_mod

            real(dp), dimension(size(phi)) :: theta
            integer :: idx

            theta = fieldline%get_theta(phi)
            do idx = 1, size(phi)
                call field%compute_B_mod(theta(idx), phi(idx), B_mod(idx))
            end do
        end subroutine B_mod_along_fieldline

    end subroutine plot_B_along_fieldline

    subroutine compare_modes(modes, labels)
        type(modes_t), dimension(:), intent(in) :: modes
        character(len=*), dimension(:), intent(in) :: labels

        type(myplot) :: plt_cos, plt_sin
        integer :: current

        call plt_cos%initialize(xlabel="mode", &
                                ylabel="cos amplitude", &
                                legend=.true.)

        do current = 1, size(modes)
            call plt_cos%add_plot(modes(current)%mode_numbers, &
                                  modes(current)%cos_coeffs, &
                                  label=labels(current), &
                                  linestyle="-o")
        end do

        call plt_cos%show()

        call plt_sin%initialize(xlabel="mode", &
                                ylabel="sin amplitude", &
                                legend=.true.)

        do current = 1, size(modes)
            call plt_sin%add_plot(modes(current)%mode_numbers, &
                                  modes(current)%sin_coeffs, &
                                  label=labels(current), &
                                  linestyle="-o")
        end do

        call plt_sin%show()

    end subroutine compare_modes

end module plot_quantities
