module shaing_callen_integration
    use constants, only: dp, pi
    use fieldline_mod, only: fieldline_t
    use integrate, only: sum_trapez_1d
    use utils, only: linspace
    implicit none
    private

    public :: get_eta_integration_grid, integrate_over_eta_grid
    public :: get_phi_integration_grid, integrate_over_phi_grid

contains

    function get_eta_integration_grid(eta_b, n_eta) result(eta_grid)
        real(dp), intent(in) :: eta_b
        integer, intent(in) :: n_eta
        real(dp), dimension(:), allocatable :: eta_grid

        real(dp), dimension(:), allocatable :: t

        allocate (eta_grid(n_eta), t(n_eta + 1))
        ! transform integral to variable t -> eta = eta_b - t**4
        ! In this transformed integral, the here considered integrands become
        ! zero at the upper limit eta = eta_b (t=0), due to the differential
        ! dt/deta = 4t**3 killing the logarithmic divergences. We therefore
        ! exlude eta = eta_b, t = 0, from the integral grid to avoid numerical
        ! instabilities due to the divergence.
        call linspace(eta_b**0.25_dp, 0.0_dp, n_eta + 1, t)
        eta_grid = eta_b - t(1:n_eta)**4.0_dp
        deallocate (t)
    end function get_eta_integration_grid

    function integrate_over_eta_grid(eta_grid, integrand) result(integral)
        real(dp), dimension(:), intent(in) :: eta_grid
        real(dp), dimension(:), intent(in) :: integrand
        real(dp) :: integral

        real(dp), parameter :: tol = 1e-15
        integer :: n_eta
        real(dp), dimension(:), allocatable :: d_eta_dt
        real(dp) :: eta_b, dt

        n_eta = size(eta_grid)

        if (n_eta < 2) then
            print *, "Error in integrate_over_eta_grid: ", &
                "need at least 2 eta points."
            error stop
        end if
        if (size(integrand) /= n_eta) then
            print *, "Error in integrate_over_eta_grid: ", &
                "size mismatch between eta_grid and integrand."
            error stop
        end if

        ! As last eta value = eta_b (1 - 1/n_eta**4) according to above grid
        eta_b = n_eta**4.0_dp/(n_eta**4.0_dp - 1.0_dp)*eta_grid(n_eta)
        if (eta_b <= maxval(eta_grid)) then
            print *, "Error in integrate_over_eta_grid: ", &
                "reconstructed eta_b not larger than all eta_grid."
            print *, "eta_b = ", eta_b, ", max(eta_grid) = ", maxval(eta_grid)
            error stop
        end if
        if (eta_b <= 0.0_dp) then
            print *, "Error in integrate_over_eta_grid: ", &
                "reconstructed eta_b must be postitiv."
            print *, "eta_b = ", eta_b
            error stop
        end if
        if (any(eta_grid < -tol*eta_b)) then
            print *, "Error in integrate_over_eta_grid: ", &
                "eta_grid must be non-negativ."
            print *, "eta_grid = ", eta_grid
            error stop
        end if

        allocate (d_eta_dt(n_eta))
        d_eta_dt = 4.0_dp*(eta_b - eta_grid)**0.75_dp
        dt = eta_b**0.25_dp/n_eta
        ! The considered integrands vanish at the upper limit (t=0) due to the
        ! differential 4t**3 (see comments in get_eta_integration_grid).
        ! As the grid is equidistant in t, the trapez sum
        !
        ! integral = (f_0/2 + f_N/2 + sum_i=1^{N-1} f_1)dt
        !
        ! has no f_N/2 term. Therefore, the right-endpoint (eta=eat_b, t=0) is
        ! excluded in eta_grid.
        integral = (sum(integrand*d_eta_dt) - 0.5_dp*integrand(1)*d_eta_dt(1))*dt
        deallocate (d_eta_dt)
    end function integrate_over_eta_grid

    function get_phi_integration_grid(fieldline, n_phi_in) result(phi_grid)
        type(fieldline_t), intent(in) :: fieldline
        real(dp), dimension(:), allocatable :: phi_grid
        integer, optional :: n_phi_in

        integer :: n_phi
        real(dp), dimension(:), allocatable :: t
        real(dp) :: delta_phi, middle_phi

        if (present(n_phi_in)) then
            n_phi = n_phi_in
        else
            n_phi = 200
        end if

        allocate (phi_grid(n_phi), t(n_phi))
        ! transform integral to variable t -> phi = middle + delta * cos(t)
        call linspace(0.0_dp, pi, n_phi, t)
        delta_phi = 0.5_dp*(fieldline%phi_max(1) - fieldline%phi_max(2))
        middle_phi = 0.5_dp*(fieldline%phi_max(1) + fieldline%phi_max(2))
        phi_grid = middle_phi + delta_phi*cos(t)
        deallocate (t)
    end function get_phi_integration_grid

    function integrate_over_phi_grid(phi_grid, integrand) result(integral)
        real(dp), dimension(:), intent(in) :: phi_grid
        real(dp), dimension(:), intent(in) :: integrand
        real(dp) :: integral

        real(dp) :: delta_phi
        integer :: n_phi
        real(dp), dimension(size(phi_grid)) :: dphi_dt, t

        n_phi = size(phi_grid)
        call linspace(0.0_dp, pi, n_phi, t)
        delta_phi = 0.5_dp*(phi_grid(1) - phi_grid(n_phi))
        dphi_dt = -delta_phi*sin(t)
        integral = sum_trapez_1d(t, integrand*dphi_dt)
    end function integrate_over_phi_grid

end module shaing_callen_integration
