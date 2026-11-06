program test_find_turning_points
    use constants, only: dp, pi
    use utils, only: linspace, not_same
    use mock_field, only: mock_field_t
    use grid_mod, only: fieldline_for_precession_t
    use precession_reference, only: find_turning_points
    use myplot_module, only: myplot

    implicit none

    real(dp), parameter :: B_0 = 1.0_dp, dB = -0.1_dp
    real(dp), parameter :: M_pol = 1.0_dp, N_tor = 1.0_dp
    type(mock_field_t) :: field
    type(fieldline_for_precession_t) :: fieldline
    real(dp), parameter :: phi_0 = 0.0_dp, theta_0 = 0.0_dp
    real(dp), parameter :: iota = 0.5_dp

    integer, parameter :: n_eta = 10
    real(dp), dimension(n_eta) :: eta
    real(dp) :: phi_turning(2), expected(2)
    real(dp) :: acos_temp

    real(dp), parameter :: reltol = 1.0e-4_dp
    real(dp) :: abstol
    integer :: current
    logical :: test_failed
    logical, parameter :: should_plot = .false.
    type(myplot) :: plt
    real(dp), dimension(100) :: phi, theta, B
    integer :: n_phi, idx

    test_failed = .false.

    call linspace(B_0 - abs(dB), B_0 + abs(dB), n_eta, eta)
    eta = 1.0_dp/eta
    eta(1) = eta(1)*(1.0_dp - reltol)
    eta(n_eta) = eta(n_eta)*(1.0_dp + reltol)

    call field%mock_field_init(M_pol, N_tor, B_0, dB)
    fieldline%phi_min = 0.0_dp
    fieldline%iota = iota
    fieldline%phi_0 = phi_0
    fieldline%theta_0 = theta_0
    fieldline%phi_max(1) = -2.0_dp*pi
    fieldline%phi_max(2) = 2.0_dp*pi
    abstol = abs(fieldline%phi_max(2) - fieldline%phi_max(1))*reltol

    do current = 1, n_eta
        phi_turning = find_turning_points(field, fieldline, eta(current), reltol)
        acos_temp = acos((1.0_dp/eta(current) - B_0)/dB)
        expected(1) = (M_pol*theta_0 - acos_temp)/(N_tor - M_pol*iota) + phi_0
        expected(2) = (M_pol*theta_0 + acos_temp)/(N_tor - M_pol*iota) + phi_0
        if (not_same(phi_turning, expected, abstol_in=abstol)) then
            print *, "test_find_turning_points failed: eta = ", eta(current)
            print *, "Expected: ", expected
            print *, "Got: ", phi_turning
            print *, "Error: ", abs(phi_turning - expected)
            print *, "Tolerance: ", abstol
            test_failed = .true.
        end if
    end do

    if (should_plot) then
        call plt%initialize(xlabel="phi", ylabel="B")
        call linspace(fieldline%phi_max(1), fieldline%phi_max(2), size(phi), phi)
        theta = fieldline%get_theta(phi)

        n_phi = size(phi)
        do idx = 1, n_phi
            call field%compute_B_mod(theta(idx), phi(idx), B(idx))
        end do
        call plt%add_plot(phi, B, label="B(phi)", linestyle="-")
        call plt%show()
    end if

    if (test_failed) error stop

end program test_find_turning_points
