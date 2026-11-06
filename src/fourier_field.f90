module fourier_field

    use constants, only: dp, pi
    use field_base, only: field_t
    use interpolate, only: SplineData2D, construct_splines_2d, evaluate_splines_2d_der

    implicit none
    private

    public :: fourier_field_t, fourier_field_init, ifft_modes_to_B

    type, extends(field_t) :: fourier_field_t
        !! Magnetic field from a flat list of Boozer Fourier modes,
        !! evaluated via 2D spline interpolation.
        !!
        !! The field strength is given by
        !! `B(theta, phi) = sum_k B_mn(k) * cos(m(k)*theta - nfp*n(k)*phi)`
        !! where `theta` and `phi` are Boozer angles and `nfp` the number of field
        !! periods.
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

    !> Initialise magnetic field from flat Fourier mode lists and
    !! evaluated via 2D spline interpolation.
    !!
    !! The field strength is given by
    !! `B(theta, phi) = sum_k B_mn(k) * cos(m(k)*theta - nfp*n(k)*phi)`
    !! where `theta` and `phi` are Boozer angles and `nfp` the number of field
    !! periods. During initialisation the series is evaluated on a
    !! `n_fft` x `n_fft` equidistant grid of the angles via fft, where
    !! `n_fft` = 2^k >= n_grid_in - 1 for smallest k. It then fits a 2D
    !! periodic quintic spline. Note again that `n` is considered normalized to
    !! the number of field periods `nfp`. If the user does not specificy
    !! `n_grid_in` explicitly, grid is chosen so Nynquist criterion is satisfied.
    subroutine fourier_field_init(field, m, n, B_mn, &
                                  B_theta_covariant, B_phi_covariant, &
                                  nfp, n_grid_in)
        use utils, only: linspace
        type(fourier_field_t), intent(inout) :: field
        integer, intent(in) :: m(:)
            !! poloidal mode numbers (flat array, length mn_max)
        integer, intent(in) :: n(:)
            !! toroidal mode numbers normalised to nfp (flat array, length mn_max)
        real(dp), intent(in) :: B_mn(:)
            !! Fourier coefficients of B in Tesla (flat array, length mn_max)
        real(dp), intent(in) :: B_theta_covariant
            !! covariant poloidal component of B in T*m
        real(dp), intent(in) :: B_phi_covariant
            !! covariant toroidal component of B in T*m
        integer, intent(in), optional :: nfp
            !! number of field periods
        integer, intent(in), optional :: n_grid_in
            !! spline grid points per angle (rounded up to 2^k + 1 for fft init)

        integer, parameter :: fft_default = 256, n_grid_max = 1025
        integer :: n_grid, n_fft
        real(dp), allocatable :: theta(:), phi(:), grid_B(:, :), fft_B(:, :)

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

        field%mn_max = max(maxval(abs(m)), maxval(abs(n)))
        if (.not. present(n_grid_in)) then
            n_fft = next_power_of_two(max(fft_default, 2*field%mn_max + 1))
            n_grid = n_fft + 1 ! closed grid for periodic spline
            if (n_grid > n_grid_max) then
                print *, "Error: Too high Fourier mode to be resolved ", &
                    "by maximum grid size!"
                print *, "m_max = ", maxval(abs(m))
                print *, "n_max = ", maxval(abs(n))
                print *, "needed n_grid = ", n_grid, " > n_grid_max = ", n_grid_max
                error stop
            end if
        elseif (n_grid_in >= 2 .and. n_grid_in <= n_grid_max) then
            n_fft = next_power_of_two(n_grid_in - 1)
            n_grid = n_fft + 1 ! closed grid for periodic spline
        else
            print *, "Error: n_grid must be between 2 and n_grid_max = ", n_grid_max
            print *, "n_grid_in = ", n_grid_in
            error stop
        end if
        field%n_grid = n_grid

        allocate (theta(n_grid), phi(n_grid))
        allocate (grid_B(n_grid, n_grid))

        call linspace(0.0_dp, 2.0_dp*pi, n_grid, theta)
        call linspace(0.0_dp, 2.0_dp*pi/field%nfp, n_grid, phi)

        allocate (fft_B(n_fft, n_fft))
        call ifft_modes_to_B(m, n, B_mn, n_fft, fft_B)
        grid_B(1:n_fft, 1:n_fft) = fft_B
        grid_B(n_grid, 1:n_fft) = grid_B(1, 1:n_fft)
        grid_B(:, n_grid) = grid_B(:, 1)

        call construct_splines_2d( &
            x_min=[0.0_dp, 0.0_dp], &
            x_max=[2.0_dp*pi, 2.0_dp*pi/field%nfp], &
            y=grid_B, &
            order=[5, 5], &
            periodic=[.true., .true.], &
            spl=field%spl)

        field%initialized = .true.

    end subroutine fourier_field_init

    function next_power_of_two(n) result(p)
        integer, intent(in) :: n
        integer :: p

        if (n < 1) then
            print *, "Error: next_power_of_two: n must be positive!"
            error stop
        end if

        p = 1
        do while (p < n)
            p = 2*p
        end do
    end function next_power_of_two

    !> In-place radix-2 FFT with a positive phase kernel and no normalisation,
    !! i.e. `z_out(p) = sum_j z_in(j) * exp(2*pi*i*(j-1)*(p-1)/n)`.
    !!
    !! `size(z)` must be a power of two. Decimation in time: the input is first
    !! put into bit-reversed order, then combined by butterflies over
    !! `log2(n)` stages, each stage reusing the partial results of the previous
    !! one.
    subroutine fft_1d(z)
        complex(dp), intent(inout) :: z(:)

        integer :: n_points, i, j, k, stage, step
        complex(dp) :: twiddle, twiddle_step, butterfly

        n_points = size(z)

        ! Bit-reversal permutation.
        j = 1
        do i = 1, n_points - 1
            if (i < j) then
                butterfly = z(j)
                z(j) = z(i)
                z(i) = butterfly
            end if
            k = n_points/2
            do
                if (k < 1) exit
                if (j <= k) exit
                j = j - k
                k = k/2
            end do
            j = j + k
        end do

        ! Butterfly stages.
        stage = 1
        do while (stage < n_points)
            step = 2*stage
            twiddle_step = cmplx(cos(pi/real(stage, dp)), &
                                 sin(pi/real(stage, dp)), dp)
            twiddle = (1.0_dp, 0.0_dp)
            do k = 1, stage
                do i = k, n_points, step
                    butterfly = twiddle*z(i + stage)
                    z(i + stage) = z(i) - butterfly
                    z(i) = z(i) + butterfly
                end do
                twiddle = twiddle*twiddle_step
            end do
            stage = step
        end do

    end subroutine fft_1d

    !> Evaluate `B(theta, phi) = sum_k B_mn(k)*cos(m(k)*theta - nfp*n(k)*phi)`
    !! on the equidistant angle grid by a single 2D inverse FFT.
    !!
    !! On the grid `theta_j = 2*pi*(j-1)/n_fft`, `phi_l = 2*pi/nfp*(l-1)/n_fft`
    !! the phase becomes `2*pi*(m*(j-1)/n_fft - n*(l-1)/n_fft)` — `nfp` cancels
    !! and every phase is a root of unity. Writing the cosine as
    !! `(exp(i*phase) + exp(-i*phase))/2` and scattering `B_mn/2` onto both
    !! `(m, -n)` and `(-m, n)` therefore turns the mode sum into an unnormalised
    !! 2D inverse DFT of a Hermitian spectrum, whose transform is real.
    subroutine ifft_modes_to_B(m, n, B_mn, n_fft, fft_B)
        integer, intent(in) :: m(:)
            !! poloidal mode numbers (flat array)
        integer, intent(in) :: n(:)
            !! toroidal mode numbers normalised to nfp (flat array)
        real(dp), intent(in) :: B_mn(:)
            !! Fourier coefficients of B in Tesla (flat array)
        integer, intent(in) :: n_fft
            !! transform length per direction, a power of two above 2*mn_max
        real(dp), intent(out) :: fft_B(:, :)
            !! field strength on the open (periodic point excluded) fft angle grid

        complex(dp), allocatable :: spectrum(:, :), work(:)
        integer :: i_mode, i_row, i_col

        allocate (spectrum(n_fft, n_fft), work(n_fft))
        spectrum = (0.0_dp, 0.0_dp)

        do i_mode = 1, size(m)
            i_row = modulo(m(i_mode), n_fft) + 1
            i_col = modulo(-n(i_mode), n_fft) + 1
            spectrum(i_row, i_col) = spectrum(i_row, i_col) &
                                     + 0.5_dp*B_mn(i_mode)
            i_row = modulo(-m(i_mode), n_fft) + 1
            i_col = modulo(n(i_mode), n_fft) + 1
            spectrum(i_row, i_col) = spectrum(i_row, i_col) &
                                     + 0.5_dp*B_mn(i_mode)
        end do

        do i_col = 1, n_fft
            call fft_1d(spectrum(:, i_col))
        end do
        do i_row = 1, n_fft
            work = spectrum(i_row, :)
            call fft_1d(work)
            spectrum(i_row, :) = work
        end do

        fft_B = real(spectrum, dp)

    end subroutine ifft_modes_to_B

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

    !> Evaluate the spline-interpolated field strength `B_mod` at Boozer angles.
    subroutine compute_B_mod(self, theta, phi, B_mod)
        class(fourier_field_t), intent(in) :: self
        real(dp), intent(in) :: theta
            !! poloidal Boozer angle in radians
        real(dp), intent(in) :: phi
            !! toroidal Boozer angle in radians
        real(dp), intent(out) :: B_mod
            !! magnetic field strength in Tesla

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
