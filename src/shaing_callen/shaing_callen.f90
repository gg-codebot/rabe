module shaing_callen_mod
    use constants, only: dp
    use fieldline_mod, only: flock_of_fieldlines_t
    use field_base, only: field_t
    use shaing_callen_integration, only: get_eta_integration_grid
    use shaing_callen_integration, only: integrate_over_eta_grid
    implicit none

contains

    !>
    !! \brief Shaing-Callen lambda_LC: computed from the trapped particle fraction
    !!
    !! \details Following the formula of Landreman and Catto (Phys. Plasmas 2012).
    !!
    !! \param[in] n_eta number of evaluation points to compute trapped fraction integral
    !! \param[in] dr_dAtheta converting lambda_LC to be used with gradients of label r
    !<
    function calc_lambda_LC(flock, field, n_eta, dr_dAtheta) result(lambda_LC)
        type(flock_of_fieldlines_t), intent(in) :: flock
        class(field_t), intent(in) :: field
        integer, intent(in) :: n_eta
        real(dp), intent(in) :: dr_dAtheta
        real(dp) :: lambda_LC

        real(dp) :: trapped_fraction, helicity_factor
        real(dp) :: B_phi_cov, B_theta_cov

        trapped_fraction = calc_trapped_fraction(flock, field, n_eta)
        call field%get_covariant_components(B_theta_cov, B_phi_cov)
        helicity_factor = (B_phi_cov*flock%M_pol + B_theta_cov*flock%N_tor) &
                          /(flock%M_pol*flock%iota - flock%N_tor)
        lambda_LC = dr_dAtheta*helicity_factor*trapped_fraction
    end function calc_lambda_LC

    function calc_trapped_fraction(flock, &
                                   field, &
                                   n_eta) result(trapped_fraction)
        type(flock_of_fieldlines_t), intent(in) :: flock
        class(field_t), intent(in) :: field
        integer, intent(in) :: n_eta
        real(dp) :: trapped_fraction

        real(dp), dimension(:), allocatable :: eta_grid
        real(dp), dimension(:), allocatable :: avg_B_squared_over_avg_lambda
        real(dp), dimension(:), allocatable :: integrand

        eta_grid = get_eta_integration_grid(flock%eta_b, n_eta)
        allocate (avg_B_squared_over_avg_lambda(n_eta))
        allocate (integrand(n_eta))
        avg_B_squared_over_avg_lambda = calc_avg_B_squared_over_avg_lambda(field, &
                                                                           flock, &
                                                                           eta_grid)
        integrand = eta_grid*avg_B_squared_over_avg_lambda
        trapped_fraction = 1.0_dp - 0.75_dp*integrate_over_eta_grid(eta_grid, &
                                                                    integrand)
        deallocate (integrand)
        deallocate (avg_B_squared_over_avg_lambda)

    end function calc_trapped_fraction

    function calc_avg_B_squared_over_avg_lambda(field, flock, eta_grid) &
        result(avg_B_squared_over_avg_lambda)
        use shaing_callen_wrappers, only: this_field
        class(field_t), intent(in) :: field
        type(flock_of_fieldlines_t), intent(in) :: flock
        real(dp), dimension(:), intent(in) :: eta_grid
        real(dp), dimension(:), allocatable :: avg_B_squared_over_avg_lambda

        integer :: this
        integer :: n_eta, n_fieldlines
        real(dp) :: avg_well_length
        real(dp), dimension(:), allocatable :: avg_lambda_over_B_squared
        real(dp), dimension(size(eta_grid)) :: fieldline_avg_loB2

        n_eta = size(eta_grid)
        n_fieldlines = size(flock%fieldlines)
        allocate (avg_lambda_over_B_squared(n_eta))
        avg_well_length = 0.0_dp
        avg_lambda_over_B_squared = 0.0_dp
        allocate (this_field, source=field)
        do this = 1, n_fieldlines
            call calc_avg_lambda_over_B_squared(flock%fieldlines(this), eta_grid, &
                                                fieldline_avg_loB2)
            avg_lambda_over_B_squared = avg_lambda_over_B_squared + fieldline_avg_loB2
            avg_well_length = avg_well_length + &
                              flock%fieldlines(this)%phi_max(2) - &
                              flock%fieldlines(this)%phi_max(1)
        end do
        deallocate (this_field)
        allocate (avg_B_squared_over_avg_lambda(n_eta))

        if (any(avg_lambda_over_B_squared <= 0.0_dp)) then
            print *, "Error in calc_avg_B_squared_over_avg_lambda: ", &
                "average lambda/B^2 is not positiv"
            print *, "avg_lambda_over_B_squared: ", avg_lambda_over_B_squared
            error stop
        end if

        avg_B_squared_over_avg_lambda = avg_well_length/ &
                                        avg_lambda_over_B_squared

        deallocate (avg_lambda_over_B_squared)

    end function calc_avg_B_squared_over_avg_lambda

    subroutine calc_avg_lambda_over_B_squared(fieldline, eta_grid, &
                                              avg_lambda_over_B_squared)
        use integrate_substituted, only: integrate_1d_substituted
        use shaing_callen_wrappers, only: wrapper_lambda_over_B_squared
        use shaing_callen_wrappers, only: this_eta, null_eta
        use shaing_callen_wrappers, only: this_fieldline, null_fieldline
        use fieldline_mod, only: fieldline_t
        type(fieldline_t), intent(in) :: fieldline
        real(dp), dimension(:), intent(in) :: eta_grid
        real(dp), dimension(:), intent(out) :: avg_lambda_over_B_squared

        integer :: this

        this_fieldline = fieldline

        do this = 1, size(eta_grid)
            this_eta = eta_grid(this)
            call integrate_1d_substituted(wrapper_lambda_over_B_squared, &
                                          fieldline%phi_max(1), &
                                          fieldline%phi_max(2), &
                                          avg_lambda_over_B_squared(this))
            this_eta = null_eta
        end do

        this_fieldline = null_fieldline

    end subroutine calc_avg_lambda_over_B_squared

    !>
    !! \brief Preliminary proxy to estimate the deviation of Lambda_LC from the
    !! true Shaing-Callen coefficient due to non-omnigeneity.
    !!
    !! \details It does not account for bootstrap resonances.
    !!
    !! \param[in] n_eta number of evaluation points to compute trapped fraction integral
    !! \param[in] dr_dAtheta converting lambda_LC to be used with gradients of label r
    !<
    function get_non_omnigenous_remainder(flock, field, n_eta, dr_dAtheta) &
        result(remainder)
        use shaing_callen_integration, only: get_eta_integration_grid
        use shaing_callen_integration, only: integrate_over_eta_grid
        type(flock_of_fieldlines_t), intent(in) :: flock
        class(field_t), intent(in) :: field
        integer, intent(in) :: n_eta
        real(dp), intent(in) :: dr_dAtheta
        real(dp) :: remainder

        real(dp), dimension(:), allocatable :: eta_grid
        real(dp), dimension(:), allocatable :: integrand
        real(dp), dimension(:), allocatable :: avg_B_squared_over_avg_lambda

        real(dp) :: B_varphi_covariant, B_vartheta_covariant

        eta_grid = get_eta_integration_grid(flock%eta_b, n_eta)
        allocate (integrand(n_eta))
        allocate (avg_B_squared_over_avg_lambda(n_eta))
        avg_B_squared_over_avg_lambda = calc_avg_B_squared_over_avg_lambda(field, &
                                                                           flock, &
                                                                           eta_grid)
        integrand = calc_avg_normalized_lambda_dphimax_dxi0(field, &
                                                            flock, &
                                                            eta_grid)
        integrand = integrand*avg_B_squared_over_avg_lambda*eta_grid
        remainder = -(calc_avg_normalized_B_squared_dphimax_dxi0(flock) - &
                      0.75_dp*integrate_over_eta_grid(eta_grid, integrand))
        deallocate (eta_grid)
        deallocate (integrand)
        deallocate (avg_B_squared_over_avg_lambda)

        call field%get_covariant_components(B_vartheta_covariant, B_varphi_covariant)
        remainder = remainder*flock%nfp/(flock%M_pol*flock%iota - flock%N_tor)
        remainder = remainder*(B_varphi_covariant + flock%iota*B_vartheta_covariant)
        remainder = remainder*dr_dAtheta

    end function get_non_omnigenous_remainder

    function calc_avg_normalized_B_squared_dphimax_dxi0(flock) &
        result(res)
        use fourier, only: real_ft
        use fieldline_labels, only: allocate_modes
        type(flock_of_fieldlines_t), intent(in) :: flock
        real(dp) :: res

        real(dp), dimension(size(flock%fieldlines)) :: w_l, phi_l, xi_0
        real(dp), dimension(size(flock%fieldlines)) :: dphimax_dxi0, well_lengths
        real(dp) :: M_pol, nfp

        M_pol = flock%M_pol
        nfp = flock%nfp
        phi_l = flock%fieldlines%phi_max(1)
        xi_0 = flock%fieldlines%xi_0
        w_l = phi_l - M_pol/nfp*xi_0
        dphimax_dxi0 = M_pol/nfp + calc_periodic_dydx(xi_0, w_l)
        well_lengths = flock%fieldlines%phi_max(2) - flock%fieldlines%phi_max(1)
        res = sum(well_lengths*dphimax_dxi0/flock%fieldlines%B_max(1)**2.0_dp)
        res = res/sum(flock%fieldlines%integral_one_over_B_squared)
    end function calc_avg_normalized_B_squared_dphimax_dxi0

    function calc_avg_normalized_lambda_dphimax_dxi0(field, flock, eta_grid) &
        result(res)
        use fourier, only: real_ft
        use fieldline_labels, only: allocate_modes
        use shaing_callen_wrappers, only: this_field
        use fieldline_integrands, only: calc_lambda_squared
        class(field_t), intent(in) :: field
        type(flock_of_fieldlines_t), intent(in) :: flock
        real(dp), dimension(:), intent(in) :: eta_grid
        real(dp), dimension(size(eta_grid)) :: res

        real(dp), dimension(size(flock%fieldlines)) :: w_l, phi_l, xi_0, B_l
        real(dp), dimension(size(flock%fieldlines)) :: dphimax_dxi0
        real(dp), dimension(size(eta_grid)) :: lambda_l
        real(dp), dimension(size(eta_grid)) :: fieldline_avg_loB2
        real(dp) :: M_pol, nfp
        integer :: this, that

        M_pol = flock%M_pol
        nfp = flock%nfp
        phi_l = flock%fieldlines%phi_max(1)
        xi_0 = flock%fieldlines%xi_0
        w_l = phi_l - M_pol/nfp*xi_0
        dphimax_dxi0 = M_pol/nfp + calc_periodic_dydx(xi_0, w_l)
        B_l = flock%fieldlines%B_max(1)
        allocate (this_field, source=field)
        res = 0.0_dp
        do this = 1, size(flock%fieldlines)
            do that = 1, size(eta_grid)
                lambda_l(that) = sqrt(calc_lambda_squared(B_l(this), eta_grid(that)))
            end do
            call calc_avg_lambda_over_B_squared(flock%fieldlines(this), eta_grid, &
                                                fieldline_avg_loB2)
            res = res + fieldline_avg_loB2*dphimax_dxi0(this)/lambda_l
        end do
        deallocate (this_field)
        res = res/sum(flock%fieldlines%integral_one_over_B_squared)
    end function calc_avg_normalized_lambda_dphimax_dxi0

    function calc_periodic_dydx(x, y) result(dydx)
        use fourier, only: check_is_equidistant, check_has_correct_endpoints
        real(dp), dimension(:), intent(in) :: x, y
        real(dp), dimension(:), allocatable :: dydx

        real(dp), dimension(:), allocatable :: y_periodic
        real(dp) :: dx

        integer :: n

        call check_is_equidistant(x)
        call check_has_correct_endpoints(x)
        n = size(x)
        if (n < 3) then
            print *, "calc_periodic_dydx needs n > 3 points!"
            print *, "provided n: ", n
            error stop
        end if

        dx = x(2) - x(1)

        allocate (y_periodic(n + 2))
        y_periodic(1) = y(n)
        y_periodic(2:n + 1) = y
        y_periodic(n + 2) = y(1)

        allocate (dydx(n))
        dydx = 0.5_dp*(y_periodic(3:n + 2) - y_periodic(1:n))/dx
        deallocate (y_periodic)

    end function calc_periodic_dydx

end module shaing_callen_mod
