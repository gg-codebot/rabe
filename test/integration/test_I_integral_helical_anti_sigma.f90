program test_I_integral_helical_anti_sigma
    use constants, only: dp, pi
    use utils, only: not_same
    use anti_sigma_field, only: anti_sigma_field_t
    use fieldline_mod, only: flock_of_fieldlines_t
    use make_fieldline, only: make_flock_of_fieldlines

    use plot_quantities, only: plot_fieldlines_over_field
    use plot_quantities, only: plot_delta_eta
    use plot_quantities, only: plot_I
    use plot_quantities, only: plot_I_factor
    use plot_quantities, only: plot_delta_A

    implicit none

    real(dp), parameter :: M_pol = 2.0_dp, N_tor = 10.0_dp
    real(dp), parameter :: B_0 = 1.0_dp, eps_0 = -0.05, eps_1 = -0.00375_dp
    real(dp), parameter :: factor = sqrt(abs(eps_0)/(1.0_dp + abs(eps_0))) &
                           /B_0**2.0_dp/N_tor*sqrt(32.0_dp)
    real(dp), parameter :: I_0_analytic = factor*(1.0_dp + abs(eps_0)*2.0_dp/3.0_dp)
    real(dp), parameter :: eps_ratio = eps_1/abs(eps_0)
    real(dp), parameter :: I_1_analytic = factor*(-0.5_dp*eps_ratio)* &
                           (1.0_dp + 6.0_dp*abs(eps_0))
    real(dp), parameter :: B_max = B_0*(1.0_dp + abs(eps_0))
    real(dp), parameter :: delta_A_1 = 0.25_dp*eps_ratio*(1.0_dp + 6.0_dp*abs(eps_0))
    type(anti_sigma_field_t) :: field
    real(dp), dimension(2), parameter :: chi_max = [+pi, -pi]

    real(dp) :: B_max_reltol, phi_error
    real(dp), parameter :: I_mean_reltol = eps_0*eps_0*3.0_dp
    real(dp), parameter :: I_amplitude_reltol = max(eps_ratio*eps_ratio*5.0_dp, &
                                                    eps_0*eps_0*10.0_dp)
    integer, parameter :: max_n_fieldlines = 20
    real(dp), dimension(2) :: phi_max

    real(dp), parameter :: iota = 0.0_dp
    real(dp), parameter :: nfp = max(1.0_dp, N_tor)
    integer :: n_fieldlines
    type(flock_of_fieldlines_t) :: flock

    real(dp) :: deviation_A, deviation_B
    real(dp) :: covariant_factor
    real(dp) :: off_factor_A, off_factor_B

    integer :: current
    real(dp), dimension(:), allocatable :: I_shifted
    real(dp) :: I_mean, I_amplitude
    logical :: test_failed

    logical, parameter :: should_plot = .false.

    call field%anti_sigma_field_init(M_pol, N_tor, B_0, eps_0, eps_1)
    call make_flock_of_fieldlines(flock, max_n_fieldlines, iota, field, &
                                  M_pol, N_tor, nfp)
    n_fieldlines = size(flock%fieldlines)
    allocate (I_shifted(n_fieldlines))

    test_failed = .false.
    do current = 1, n_fieldlines
        phi_max = (M_pol*flock%fieldlines(current)%theta_0 - chi_max)/N_tor
        if (not_same(flock%fieldlines(current)%phi_max(1), &
                     phi_max(1), &
                     abstol_in=2.0_dp*flock%fieldlines(current)%phi_max_error(1))) then
            print *, "-------------------------------------------------------------"
            print *, "test_I_integral_helical_anti_sigma failed: left phi_max"
            print *, "found: ", flock%fieldlines(current)%phi_max(1)
            print *, "analytic: ", phi_max(1)
            test_failed = .true.
        end if
        if (not_same(flock%fieldlines(current)%phi_max(2), &
                     phi_max(2), &
                     abstol_in=2.0_dp*flock%fieldlines(current)%phi_max_error(2))) then
            print *, "-------------------------------------------------------------"
            print *, "test_I_integral_helical_anti_sigma failed: right phi_max"
            print *, "found: ", flock%fieldlines(current)%phi_max(2)
            print *, "analytic: ", phi_max(2)
            test_failed = .true.
        end if

    end do

    phi_error = 2.0_dp*max(maxval(flock%fieldlines%phi_max_error(1)), &
                           maxval(flock%fieldlines%phi_max_error(2)))

    B_max_reltol = phi_error**2.0_dp*N_tor**2.0
    if (not_same(1.0_dp/flock%eta_b, B_max, reltol_in=B_max_reltol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_I_integral_helical_anti_sigma failed: B_max"
        print *, "found: ", 1.0_dp/flock%eta_b
        print *, "analytic: ", B_max
        print *, "relative error: ", flock%eta_b*B_max - 1.0_dp
        test_failed = .true.
    end if

    I_mean = sum(flock%fieldlines%integral_lambda_b_over_B_squared)/n_fieldlines
    if (not_same(I_mean, I_0_analytic, reltol_in=I_mean_reltol)) then
        print *, "-------------------------------------------------------------"
        print *, "test_I_integral_helical_anti_sigma failed: I_mean"
        print *, "found: ", I_mean
        print *, "analytic: ", I_0_analytic
        print *, "relative error: ", I_mean/I_0_analytic - 1.0_dp
        test_failed = .true.
    end if

    I_shifted = flock%fieldlines%integral_lambda_b_over_B_squared - I_mean
    I_amplitude = 0.5_dp*(maxval(I_shifted) - minval(I_shifted))

    if (not_same(I_amplitude, I_1_analytic, &
                 reltol_in=I_amplitude_reltol, abstol_in=0.0_dp)) then
        print *, "-------------------------------------------------------------"
        print *, "test_I_integral_helical_anti_sigma failed: I_amplitude"
        print *, "found: ", I_amplitude
        print *, "analytic: ", I_1_analytic
        print *, "relative error: ", I_amplitude/I_1_analytic - 1.0_dp
        test_failed = .true.
    end if

    if (test_failed) error stop

    if (should_plot) then
        call plot_fieldlines_over_field(flock%fieldlines, field)
        call plot_delta_eta(flock%fieldlines, flock%iota_p)
        call plot_I_factor(flock%fieldlines, eps_0, eps_1)
        call plot_I(flock%fieldlines, I_0_analytic, I_1_analytic, eps_0, eps_1)
        call plot_delta_A(flock%fieldlines, delta_A_1)
    end if

end program test_I_integral_helical_anti_sigma
