module fourier_field

    use constants, only: dp, pi
    use field_base, only: field_t
    use interpolate, only: SplineData2D, construct_splines_2d, evaluate_splines_2d_der

    implicit none

    type, extends(field_t) :: fourier_field_t
        !>
        !! \brief Magnetic field from a flat list of Boozer Fourier modes,
        !! evaluated via 2D spline interpolation.
        !!
        !! \details The field strength is given by
        !! `B(theta, phi) = sum_k B_mn(k) * cos(m(k)*theta - nfp*n(k)*phi)`
        !! where `theta` and `phi` are Boozer angles and `nfp` the number of field
        !! periods.
        !<
        logical :: initialized = .false.
        type(SplineData2D) :: spl
        real(dp) :: nfp
        real(dp) :: B_theta_covariant, B_phi_covariant
        integer :: n_grid, mn_max
    contains
        procedure :: compute_B_sqrtg_dB_dx
        procedure :: compute_B_and_dB_dx
        procedure :: compute_B_mod
        procedure :: compute_nabla_s
        procedure :: rel_accuracy_B
        procedure :: get_covariant_components
    end type fourier_field_t

contains

    !>
    !! \brief Initialise magnetic field from flat Fourier mode lists and
    !! evaluated via 2D spline interpolation.
    !!
    !! \details The field strength is given by
    !! `B(theta, phi) = sum_k B_mn(k) * cos(m(k)*theta - nfp*n(k)*phi)`
    !! where `theta` and `phi` are Boozer angles and `nfp` the number of field
    !! periods. During initialisation the series is evaluated on a
    !! `n_grid` x `n_grid` equidistant grid of the angles, then fits a 2D
    !! periodic quintic spline. Note again that `n` is considered normalized to
    !! the number of field periods `nfp`.
    !!
    !! \param[in] m poloidal mode numbers (flat array, length mn_max)
    !! \param[in] n toroidal mode numbers normalised to nfp (flat array, length mn_max)
    !! \param[in] B_mn Fourier coefficients of B in Tesla (flat array, length mn_max)
    !! \param[in] B_theta_covariant covariant poloidal component of B in T*m
    !! \param[in] B_phi_covariant covariant toroidal component of B in T*m
    !! \param[in] nfp number of field periods
    !! \param[in] n_grid spline grid points per angle direction
    !<
    subroutine fourier_field_init(field, m, n, B_mn, &
                                  B_theta_covariant, B_phi_covariant, &
                                  nfp, n_grid)
        use utils, only: linspace
        type(fourier_field_t), intent(inout) :: field
        integer, intent(in) :: m(:), n(:)
        real(dp), intent(in) :: B_mn(:)
        real(dp), intent(in) :: B_theta_covariant, B_phi_covariant
        integer, intent(in), optional :: nfp
        integer, intent(in), optional :: n_grid

        integer, parameter :: n_grid_default = 200

        integer :: n_theta, n_phi, i_theta, i_phi
        real(dp), allocatable :: theta(:), phi(:), grid_B(:, :)

        field%B_theta_covariant = B_theta_covariant
        field%B_phi_covariant = B_phi_covariant

        if (.not. present(nfp)) then
            field%nfp = 1.0_dp
        elseif (nfp > 0) then
            field%nfp = real(nfp, dp)
        else
            print *, "Error: nfp must be positive!"
            error stop
        end if

        if (.not. present(n_grid)) then
            n_theta = n_grid_default
            n_phi = n_grid_default
        elseif (n_grid >= 2) then
            n_theta = n_grid
            n_phi = n_grid
        else
            print *, "Error: n_grid must be at least 2!"
            error stop
        end if
        field%n_grid = n_theta
        field%mn_max = max(maxval(abs(m)), maxval(abs(n)))

        allocate (theta(n_theta), phi(n_phi))
        allocate (grid_B(n_theta, n_phi))

        call linspace(0.0_dp, 2.0_dp*pi, n_theta, theta)
        call linspace(0.0_dp, 2.0_dp*pi/field%nfp, n_phi, phi)

        do i_phi = 1, n_phi
            do i_theta = 1, n_theta
                grid_B(i_theta, i_phi) = sum(B_mn*cos(real(m, dp)*theta(i_theta) &
                                                    - real(n, dp)*field%nfp*phi(i_phi)))
            end do
        end do

        call construct_splines_2d( &
            x_min=[0.0_dp, 0.0_dp], &
            x_max=[2.0_dp*pi, 2.0_dp*pi/field%nfp], &
            y=grid_B, &
            order=[5, 5], &
            periodic=[.true., .true.], &
            spl=field%spl)

        field%initialized = .true.

    end subroutine fourier_field_init

    subroutine compute_B_sqrtg_dB_dx(self, theta, phi, B_mod, sqrtg, dB_dx)
        class(fourier_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod, sqrtg, dB_dx(3)

        if (.not. self%initialized) &
            error stop "compute_B_sqrtg_dB_dx: fourier_field_t not initialized"
        print *, "fourier_field_t does not provide sqrtg. "// &
            "Use compute_B_mod or compute_B_and_dB_dx instead!"
        error stop

    end subroutine compute_B_sqrtg_dB_dx

    subroutine compute_B_and_dB_dx(self, theta, phi, B_mod, dB_dx)
        class(fourier_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod, dB_dx(3)

        real(dp) :: x(2), dy(2)

        if (.not. self%initialized) &
            error stop "compute_B_and_dB_dx: fourier_field_t not initialized"
        x(1) = theta
        x(2) = phi
        call evaluate_splines_2d_der(self%spl, x, B_mod, dy)
        dB_dx(1) = 0.0_dp
        dB_dx(2) = dy(1)
        dB_dx(3) = dy(2)

    end subroutine compute_B_and_dB_dx

    !>
    !! \brief Evaluate the spline-interpolated field strength `B_mod` at Boozer angles.
    !!
    !! \param[in] theta poloidal Boozer angle in radians
    !! \param[in] phi toroidal Boozer angle in radians
    !! \param[out] B_mod magnetic field strength in Tesla
    !<
    subroutine compute_B_mod(self, theta, phi, B_mod)
        class(fourier_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod

        real(dp) :: dummy(3)

        call self%compute_B_and_dB_dx(theta, phi, B_mod, dummy)

    end subroutine compute_B_mod

    subroutine compute_nabla_s(self, theta, phi, nabla_s)
        class(fourier_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: nabla_s

        if (.not. self%initialized) &
            error stop "compute_nabla_s: fourier_field_t not initialized"
        nabla_s = 1.0_dp

    end subroutine compute_nabla_s

    real(dp) function rel_accuracy_B(self)
        class(fourier_field_t), intent(in) :: self

        rel_accuracy_B = (real(self%mn_max, dp)/real(self%n_grid, dp))**6
    end function rel_accuracy_B

    subroutine get_covariant_components(self, B_theta_covariant, B_phi_covariant)
        class(fourier_field_t), intent(in) :: self
        real(dp), intent(out) :: B_theta_covariant, B_phi_covariant

        if (.not. self%initialized) &
            error stop "get_covariant_components: fourier_field_t not initialized"
        B_theta_covariant = self%B_theta_covariant
        B_phi_covariant = self%B_phi_covariant

    end subroutine get_covariant_components

end module fourier_field
