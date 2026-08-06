module grid_mod
    use constants, only: dp, pi
    use field_base, only: field_3D_t
    use fieldline_mod, only: fieldline_t
    use interpolate, only: SplineData1D, construct_splines_1d, evaluate_splines_1d_many
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    implicit none

    type :: integration_grid_t
        integer :: n
        real(dp), dimension(:), allocatable :: eta
        real(dp), dimension(:), allocatable :: t
        real(dp), dimension(:), allocatable :: normalized_bounce_time
        real(dp), dimension(:), allocatable :: bounce_time_weighted
        real(dp), dimension(:), allocatable :: I_j
        real(dp), dimension(:), allocatable :: normalized_radial_drift
        real(dp), dimension(:), allocatable :: radial_drift_weighted
        real(dp), dimension(:), allocatable :: poloidal_drift
        real(dp), dimension(:), allocatable :: poloidal_drift_weighted
        type(SplineData1D) :: I_j_spline
        type(SplineData1D) :: bounce_time_weighted_spline
        type(SplineData1D) :: radial_drift_weighted_spline
        type(SplineData1D) :: poloidal_drift_weighted_spline
    end type integration_grid_t

    type, extends(fieldline_t) :: fieldline_for_precession_t
        real(dp) :: M_pol
        real(dp) :: N_tor
        real(dp) :: nfp
        real(dp) :: phi_min
        real(dp) :: B_min
        type(integration_grid_t) :: grid
    end type fieldline_for_precession_t

contains

    subroutine set_integration_grids(eta_t, eta_c, grid, n_eta)
        use utils, only: linspace
        use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
        real(dp), intent(in) :: eta_t, eta_c
        type(integration_grid_t), intent(out) :: grid
        integer, intent(in), optional :: n_eta

        real(dp) :: eta_mid
        real(dp) :: t_start, t_end
        integer :: n
        real(dp), dimension(:), allocatable :: t

        if (ieee_is_nan(eta_t)) then
            print *, "Error: eta_t is NaN."
            error stop
        end if
        if (ieee_is_nan(eta_c)) then
            print *, "Error: eta_c is NaN."
            error stop
        end if
        if (eta_t <= eta_c) then
            print *, "Error: eta_t must be greater than eta_c."
            error stop
        end if

        if (present(n_eta)) then
            if (n_eta < 2) then
                print *, "Error: n_eta must be at least 2."
                error stop
            end if
            n = n_eta
        else
            n = 200
        end if

        if (allocated(grid%eta)) deallocate (grid%eta)
        if (allocated(grid%t)) deallocate (grid%t)
        if (allocated(t)) deallocate (t)

        allocate (grid%eta(n))
        allocate (grid%t(n))
        allocate (t(n))
        t_start = (eta_t - eta_c)**0.5_dp
        t_end = 0.0_dp
        call linspace(t_start, t_end, n, t)
        grid%eta = eta_c + t**2.0_dp
        grid%t = t
        grid%n = n

        if (allocated(grid%bounce_time_weighted)) then
            deallocate (grid%bounce_time_weighted)
        end if
        if (allocated(grid%normalized_bounce_time)) then
            deallocate (grid%normalized_bounce_time)
        end if
        if (allocated(grid%I_j)) then
            deallocate (grid%I_j)
        end if
        if (allocated(grid%normalized_radial_drift)) then
            deallocate (grid%normalized_radial_drift)
        end if
        if (allocated(grid%radial_drift_weighted)) then
            deallocate (grid%radial_drift_weighted)
        end if
        if (allocated(grid%poloidal_drift)) then
            deallocate (grid%poloidal_drift)
        end if
        if (allocated(grid%poloidal_drift_weighted)) then
            deallocate (grid%poloidal_drift_weighted)
        end if

        allocate (grid%bounce_time_weighted(n))
        allocate (grid%normalized_bounce_time(n))
        allocate (grid%I_j(n))
        allocate (grid%normalized_radial_drift(n))
        allocate (grid%radial_drift_weighted(n))
        allocate (grid%poloidal_drift(n))
        allocate (grid%poloidal_drift_weighted(n))

        grid%bounce_time_weighted = ieee_value(1.0_dp, ieee_quiet_nan)
        grid%normalized_bounce_time = ieee_value(1.0_dp, ieee_quiet_nan)
        grid%I_j = ieee_value(1.0_dp, ieee_quiet_nan)
        grid%normalized_radial_drift = ieee_value(1.0_dp, ieee_quiet_nan)
        grid%radial_drift_weighted = ieee_value(1.0_dp, ieee_quiet_nan)
        grid%poloidal_drift = ieee_value(1.0_dp, ieee_quiet_nan)
        grid%poloidal_drift_weighted = ieee_value(1.0_dp, ieee_quiet_nan)

        if (any(ieee_is_nan(grid%eta))) then
            print *, "Error: NaN values found in eta grid."
            error stop
        end if
        if (any(ieee_is_nan(grid%t))) then
            print *, "Error: NaN values found in t grid."
            error stop
        end if

    end subroutine set_integration_grids

    subroutine compute_bounce_integrals(field, fieldline, s_tor, grid, orbit_history)
        use field_instance, only: initialize_field_instance
        use params, only: params_init, ntau, ntimstep
        use bounce, only: trace_orbit_till_bounce
        use fieldline_integrands, only: calc_lambda_squared
        use constants, only: machine_eps
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        class(field_3D_t), intent(in) :: field
        type(fieldline_for_precession_t), intent(in) :: fieldline
        real(dp), intent(in) :: s_tor
        type(integration_grid_t), intent(inout) :: grid
     real(dp), dimension(:, :, :), allocatable, intent(inout), optional :: orbit_history

        integer, parameter :: n_coef = 2
        integer :: idx, n_grid

        real(dp) :: theta_min, phi_min

        integer, parameter :: n_dim = 7
        real(dp), parameter :: cm2m = 1e-2_dp, gauss2tesla = 1e-4_dp
        real(dp), parameter :: lamor_radius = 1e-8_dp !m*Tesla
        real(dp) :: well_depth
        real(dp) :: deep_trapped_bounce_time, time_step
        real(dp), dimension(n_dim) :: z_template, z_start, z_end
        real(dp) :: lambda_squared
        real(dp) :: eta, t_weight
        real(dp) :: bounce_time_times_v_thermal
        real(dp) :: I_j
        real(dp) :: h_ctrvr(3), x(3)
        real(dp) :: dummy(11)
        integer :: max_trace_steps

        theta_min = fieldline%get_theta(fieldline%phi_min)
        phi_min = fieldline%phi_min
        x = [s_tor, theta_min, phi_min]
        call field%evaluate(x, dummy(1), dummy(2), dummy(3:5), dummy(6:8), &
                            h_ctrvr, dummy(9:11))

        well_depth = 1.0_dp - fieldline%B_min*fieldline%eta_b
        deep_trapped_bounce_time = 4.0_dp*pi/abs(h_ctrvr(2))/sqrt(well_depth) &
                                  /abs(fieldline%M_pol*fieldline%iota - fieldline%N_tor)
        time_step = deep_trapped_bounce_time/50.0_dp
        call params_init(fieldline%nfp, &
                         time_step, &
                         rlarm_in=lamor_radius)
        max_trace_steps = ntau*(ntimstep - 1)
        if (present(orbit_history)) then
            if (.not. allocated(orbit_history)) then
                allocate (orbit_history(max_trace_steps, n_dim, grid%n))
            end if
            orbit_history = ieee_value(1.0_dp, ieee_quiet_nan)
        end if
        z_template(1) = s_tor
        z_template(2) = theta_min
        z_template(3) = phi_min
        z_template(4) = 1.0_dp
        z_template(5) = 0.0_dp
        z_template(6) = 0.0_dp
        z_template(7) = 0.0_dp

        do idx = 1, grid%n
            t_weight = 2.0_dp*grid%t(idx)
            if (abs(t_weight) < machine_eps) then
                grid%bounce_time_weighted(idx) = 0.0_dp
                grid%normalized_bounce_time(idx) = ieee_value(1.0_dp, ieee_quiet_nan)
                grid%I_j(idx) = ieee_value(1.0_dp, ieee_quiet_nan)
                grid%normalized_radial_drift(idx) = ieee_value(1.0_dp, ieee_quiet_nan)
                grid%radial_drift_weighted(idx) = 0.0_dp
                grid%poloidal_drift(idx) = ieee_value(1.0_dp, ieee_quiet_nan)
                grid%poloidal_drift_weighted(idx) = 0.0_dp
                cycle
            end if

            !> initialize orbit with certain eta
            z_start = z_template
            eta = grid%eta(idx)
            lambda_squared = calc_lambda_squared(fieldline%B_min, eta)
            z_start(5) = sqrt(lambda_squared)
            if (ieee_is_nan(z_start(5))) then
                print *, "Error: Lambda is NaN!"
                print *, "B_min:", fieldline%B_min, "eta:", eta
                print *, "Calculated lambda_squared:", lambda_squared
                error stop
            end if

            !> initial trace to a turning point
            call trace_orbit_till_bounce(z_start, z_end)

            !> then reset time and integral and trace to until back at bounce
            !> i.e. one full bounce period
            z_start = z_end
            z_start(6:7) = 0.0_dp
            if (present(orbit_history)) then
                call trace_orbit_till_bounce(z_start, z_end, orbit_history(:, :, idx))
            else
                call trace_orbit_till_bounce(z_start, z_end)
            end if
            I_j = 0.5_dp*z_end(6)*cm2m/gauss2tesla ! I_j is only over half a bounce period
            bounce_time_times_v_thermal = z_end(7)*cm2m
            grid%normalized_bounce_time(idx) = bounce_time_times_v_thermal
            grid%bounce_time_weighted(idx) = t_weight*grid%normalized_bounce_time(idx)
            grid%I_j(idx) = I_j
            grid%normalized_radial_drift(idx) = (z_end(1) - z_start(1))/lamor_radius
            grid%radial_drift_weighted(idx) = t_weight*grid%normalized_radial_drift(idx)
            grid%poloidal_drift(idx) = (z_end(2) - z_start(2))/lamor_radius
            grid%poloidal_drift_weighted(idx) = t_weight*grid%poloidal_drift(idx)

            ! write (out_unit, '(I0,1X,8(ES22.14E3,1X))') &
            !     idx, grid%t(idx), eta, &
            !     bounce_time_times_v_thermal, I_j, &
            !     (z_end(1) - z_start(1))/lamor_scan(i_lamor), &
            !     t_weight*((z_end(1) - z_start(1))/lamor_scan(i_lamor)), &
            !     (z_end(2) - z_start(2))/lamor_scan(i_lamor), &
            !     t_weight*((z_end(2) - z_start(2))/lamor_scan(i_lamor))

            print *, "idx:", idx, "out of", grid%n
        end do

    end subroutine compute_bounce_integrals

    subroutine set_splines(grid)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
        type(integration_grid_t), intent(inout) :: grid
        logical, parameter :: periodic = .false.
        integer, parameter :: I_j_order = 3
        integer, parameter :: bounce_order = 3

        integer :: n
        integer :: start

        if (ieee_is_nan(grid%I_j(1))) then
            start = 2
        else
            start = 1
        end if
        if (ieee_is_nan(grid%I_j(grid%n))) then
            n = grid%n - 1
        else
            n = grid%n
        end if

        call construct_splines_1d(grid%t(start), &
                                  grid%t(n), &
                                  grid%I_j(start:n), &
                                  I_j_order, &
                                  periodic, &
                                  grid%I_j_spline)

        start = 1
        n = grid%n
        call construct_splines_1d(grid%t(start), &
                                  grid%t(n), &
                                  grid%bounce_time_weighted(start:n), &
                                  bounce_order, &
                                  periodic, &
                                  grid%bounce_time_weighted_spline)

        call construct_splines_1d(grid%t(start), &
                                  grid%t(n), &
                                  grid%radial_drift_weighted(start:n), &
                                  bounce_order, &
                                  periodic, &
                                  grid%radial_drift_weighted_spline)

        call construct_splines_1d(grid%t(start), &
                                  grid%t(n), &
                                  grid%poloidal_drift_weighted(start:n), &
                                  bounce_order, &
                                  periodic, &
                                  grid%poloidal_drift_weighted_spline)
    end subroutine set_splines

    subroutine evaluate_grid_splines(grid, t_eval, I_j_eval, bounce_coef_eval, &
                                     radial_drift_eval, poloidal_drift_eval)
        type(integration_grid_t), intent(in) :: grid
        real(dp), intent(in) :: t_eval(:)
        real(dp), intent(out) :: I_j_eval(:)
        real(dp), intent(out) :: bounce_coef_eval(:)
        real(dp), intent(out), optional :: radial_drift_eval(:)
        real(dp), intent(out), optional :: poloidal_drift_eval(:)

        call evaluate_splines_1d_many(grid%I_j_spline, t_eval, I_j_eval)
        call evaluate_splines_1d_many(grid%bounce_time_weighted_spline, t_eval, bounce_coef_eval)
        if (present(radial_drift_eval)) then
            call evaluate_splines_1d_many(grid%radial_drift_weighted_spline, t_eval, radial_drift_eval)
        end if
        if (present(poloidal_drift_eval)) then
            call evaluate_splines_1d_many(grid%poloidal_drift_weighted_spline, t_eval, poloidal_drift_eval)
        end if
    end subroutine evaluate_grid_splines
end module grid_mod
