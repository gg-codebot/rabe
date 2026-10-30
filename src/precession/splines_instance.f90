
module splines_instance
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use constants, only: dp, machine_eps
    use interpolate, only: evaluate_splines_1d
    use interpolate, only: SplineData1D
    use interpolate, only: construct_splines_1d

    implicit none

    private
    type(SplineData1D) :: I_j_spl
    type(SplineData1D) :: poloidal_drift_spl
    type(SplineData1D) :: radial_drift_mode_spl
    type(SplineData1D) :: eta_spl
    real(dp) :: prefactor

    logical :: splines_initialized = .false.
    logical, parameter :: periodic = .false.
    integer, parameter :: order = 3 ! must stay 3

    type(SplineData1D) :: startup_spl
    complex(dp), parameter :: i = (0.0_dp, 1.0_dp)
    complex(dp), dimension(order) :: y_hat

    public :: get_flux_mode
    public :: get_radial_drift_mode
    public :: get_poloidal_drift
    public :: get_I_j
    public :: initialize_splines
    public :: initialize_radial_drift_spline
    public :: initialize_prefactor
    public :: initialize_startup

contains
    subroutine get_flux_mode(t_start, t_end, flux_mode)
        use odeint_allroutines_sub, only: odeint_allroutines
        real(dp), intent(in) :: t_start, t_end
        real(dp), intent(out) :: flux_mode

        real(dp), parameter :: relerr = 1.0e-8_dp
        integer, parameter :: ndim = 6
        real(dp) :: y(ndim), y_startup(ndim)
        real(dp), parameter :: rel_start_up = 1.0e-4_dp
        complex(dp) :: scaler
        real(dp) :: eps_t
        real(dp) :: I_j, eta
        if (.not. splines_initialized) then
          print *, "Error: get_flux_mode called before splines instance is initialized."
            error stop
        end if

        eps_t = rel_start_up*(t_end - t_start)
        call get_startup(y_startup, t_start, eps_t)
        y = y_startup
        call odeint_allroutines(y, ndim, t_start + eps_t, t_end, relerr, rhs)

        scaler = -i/complex(y(3), -y(1))
        flux_mode = real(i*scaler*complex(y(6), -y(5)))
    end subroutine get_flux_mode

    subroutine get_startup(y_startup, t_start, eps_t)
        real(dp), dimension(:), intent(out) :: y_startup
        real(dp), intent(in) :: t_start
        real(dp), intent(in) :: eps_t

        real(dp) :: I_j, eta
        complex(dp) :: complex_y(2)

        call get_I_j(t_start + eps_t, I_j)
        call get_eta(t_start + eps_t, eta)
        complex_y(1) = (1.0_dp, 0.0_dp) + &
            y_hat(1)*eps_t + &
            y_hat(2)*eps_t**2 + &
            y_hat(3)*eps_t**3
        complex_y(2) = y_hat(1) + &
            2.0_dp*y_hat(2)*eps_t + &
            3.0_dp*y_hat(3)*eps_t**2

        y_startup(1) = -aimag(complex_y(1))
        y_startup(2) = -aimag(complex_y(2))*eta*I_j/(2.0_dp*(t_start + eps_t))
        y_startup(3) = real(complex_y(1))
        y_startup(4) = real(complex_y(2))*eta*I_j/(2.0_dp*(t_start + eps_t))
        y_startup(5) = 0.0_dp
        y_startup(6) = 0.0_dp
    end subroutine get_startup

    subroutine rhs(t, y, dydt)
        real(dp), intent(in) :: t
        real(dp), dimension(:), intent(in) :: y
        real(dp), dimension(:), intent(out) :: dydt

        real(dp) :: radial_drift_mode, poloidal_drift, I_j
        real(dp) :: eta

        call get_radial_drift_mode(t, radial_drift_mode)
        call get_poloidal_drift(t, poloidal_drift)
        call get_I_j(t, I_j)
        call get_eta(t, eta)

        dydt(1) = 2.0_dp*t/(I_j*eta)*y(2)
        dydt(2) = -poloidal_drift*prefactor*y(3)
        dydt(3) = 2.0_dp*t/(I_j*eta)*y(4)
        dydt(4) = poloidal_drift*prefactor*y(1)
        dydt(5) = radial_drift_mode*y(1)
        dydt(6) = radial_drift_mode*y(3)

    end subroutine rhs

    subroutine get_radial_drift_mode(t, radial_drift)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: radial_drift

        if (.not. splines_initialized) then
     error stop "Error: get_radial_drift called before splines instance is initialized."
        end if
        call evaluate_splines_1d(radial_drift_mode_spl, t, radial_drift)
    end subroutine get_radial_drift_mode

    subroutine get_poloidal_drift(t, magnetic_drift)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: magnetic_drift

        if (.not. splines_initialized) then
   error stop "Error: get_magnetic_drift called before splines instance is initialized."
        end if
        call evaluate_splines_1d(poloidal_drift_spl, t, magnetic_drift)
    end subroutine get_poloidal_drift

    subroutine get_I_j(t, I_j)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: I_j

        if (.not. splines_initialized) then
            error stop "Error: get_I_j called before splines instance is initialized."
        end if
        call evaluate_splines_1d(I_j_spl, t, I_j)
    end subroutine get_I_j

    subroutine get_eta(t, eta)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: eta

        if (.not. splines_initialized) then
            print *, "Error: get_eta called before splines instance is initialized."
            error stop
        end if
        call evaluate_splines_1d(eta_spl, t, eta)
    end subroutine get_eta

    subroutine initialize_splines(t, eta, I_j, ploidal_drift)
        real(dp), dimension(:), intent(in) :: t
        real(dp), dimension(:), intent(in) :: eta
        real(dp), dimension(:), intent(in) :: I_j
        real(dp), dimension(:), intent(in) :: ploidal_drift

        integer :: start, n

        n = size(t)
        if (ieee_is_nan(I_j(1))) then
            start = 2
        else
            start = 1
        end if
        if (ieee_is_nan(I_j(n))) then
            n = size(t) - 1
        else
            n = size(t)
        end if
        if (any(ieee_is_nan(I_j(start:n)))) then
            print *, "I_j: ", I_j(start:n)
            error stop "Error: initialize_splines called with NaN in I_j."
        end if
        call construct_splines_1d(t(start), t(n), I_j(start:n), &
                                  order, periodic, I_j_spl)

        start = 1
        n = size(t)
        if (any(ieee_is_nan(ploidal_drift(start:n)))) then
            print *, "poloidal_drift: ", ploidal_drift(start:n)
            error stop "Error: initialize_splines called with NaN in poloidal_drift."
        end if
        call construct_splines_1d(t(start), t(n), ploidal_drift(start:n), &
                                  order, periodic, poloidal_drift_spl)
        if (any(ieee_is_nan(eta(start:n)))) then
            print *, "eta: ", eta(start:n)
            error stop "Error: initialize_splines called with NaN in eta."
        end if
        call construct_splines_1d(t(start), t(n), eta(start:n), &
                                  order, periodic, eta_spl)

        splines_initialized = .true.
    end subroutine initialize_splines

    subroutine initialize_startup(t, eta, I_j, mode_factor)
        real(dp), dimension(:), intent(in) :: t
        real(dp), dimension(:), intent(in) :: eta
        real(dp), dimension(:), intent(in) :: I_j
        real(dp), intent(in) :: mode_factor

        integer :: start, n
        real(dp) :: t_mid
        real(dp), dimension(:), allocatable :: startup
        real(dp), dimension(0:order) :: a, b

        if (any(ieee_is_nan(t))) then
            error stop "Error: initialize_startup called with NaN in t."
        end if
        if (any(ieee_is_nan(eta))) then
            error stop "Error: initialize_startup called with NaN in eta."
        end if
        if (ieee_is_nan(mode_factor)) then
            error stop "Error: initialize_startup called with NaN in mode_factor."
        end if
        if (abs(mode_factor) < machine_eps) then
            print *, "Error: initialize_startup called with mode_factor too ", &
                "close to zero."
            error stop
        end if

        t_mid = 0.5_dp*(t(1) + t(size(t)))
        n = count(t < t_mid)
        start = 1
        if (any(ieee_is_nan(I_j(start:n)))) then
            error stop "Error: initialize_startup called with NaN in I_j."
        end if
        if (any(abs(t(start:n)) < machine_eps)) then
            print *, "t = ", t(start:n)
            print *, "Error: initialize_startup called with t values too close to zero."
            error stop
        end if

        allocate (startup(n))
        startup = 0.5_dp*I_j(start:n)*eta(start:n)/t(start:n)/mode_factor
        call construct_splines_1d(t(start), t(n), startup, &
                                  order, periodic, startup_spl)
        deallocate (startup)
        a = poloidal_drift_spl%coeff(:, 1)
        b = startup_spl%coeff(:, 1)
        y_hat(1) = i*a(0)/b(1)
        y_hat(2) = ((i*a(0) - 2.0_dp*b(2))*y_hat(1) + i*a(1))/(4.0_dp*b(1))
        y_hat(3) = ((i*a(0) - 6.0_dp*b(2))*y_hat(2) + &
                    (i*a(1) - 3.0_dp*b(2))*y_hat(1) + i*a(2))/ &
                   (9.0_dp*b(1))
    end subroutine initialize_startup

    subroutine initialize_radial_drift_spline(t, radial_drift_mode)
        real(dp), dimension(:), intent(in) :: t
        real(dp), dimension(:), intent(in) :: radial_drift_mode

        integer :: start, n

        n = size(t)
        start = 1
        if (any(ieee_is_nan(radial_drift_mode(start:n)))) then
            print *, "radial_drift_mode: ", radial_drift_mode(start:n)
            print *, "Error: initialize_radial_drift_spline called with NaN in ", &
                "radial_drift_mode."
            error stop
        end if
        call construct_splines_1d(t(start), t(n), radial_drift_mode(start:n), &
                                  order, periodic, radial_drift_mode_spl)
    end subroutine initialize_radial_drift_spline

    subroutine initialize_prefactor(prefactor_in)
        real(dp), intent(in) :: prefactor_in
        if (ieee_is_nan(prefactor_in)) then
            print *, "prefactor_in: ", prefactor_in
            error stop "Error: initialize_prefactor called with NaN in prefactor_in."
        end if
        prefactor = prefactor_in
    end subroutine initialize_prefactor
end module splines_instance
