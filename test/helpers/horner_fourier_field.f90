module horner_fourier_field

    use constants, only: dp
    use field_base, only: field_t

    implicit none

    type, extends(field_t) :: horner_fourier_field_t
        integer :: number_of_modes
        integer, allocatable :: m(:), n(:)
        real(dp), allocatable :: B_mn(:)
        integer :: nfp, n_min
    contains
        procedure :: horner_fourier_field_init
        procedure :: compute_B_sqrtg_dB_dx
        procedure :: compute_B_and_dB_dx
        procedure :: compute_B_mod
        procedure :: compute_nabla_s
        procedure :: rel_accuracy_B
        procedure :: get_covariant_components
    end type horner_fourier_field_t

contains

    ! Initialize from 2D coefficient array.
    ! B = sum_{m,n} B_mn(m_idx,n_idx) * cos(m*theta - n*phi)
    subroutine horner_fourier_field_init(self, m_arr, n_arr, B_mn, nfp)
        class(horner_fourier_field_t), intent(out) :: self
        integer, intent(in) :: m_arr(:), n_arr(:)
        real(dp), intent(in) :: B_mn(:)
        integer, intent(in), optional :: nfp

        integer :: nm, nn, im, in_, k

        integer, allocatable :: idx(:)
        integer :: i, j, n
        integer :: temp_idx

        self%number_of_modes = size(m_arr)
        if (present(nfp)) then
            self%nfp = nfp
        else
            self%nfp = 1
        end if

        allocate (self%m(self%number_of_modes))
        allocate (self%n(self%number_of_modes))
        allocate (self%B_mn(self%number_of_modes))

        allocate (idx(self%number_of_modes))
        do i = 1, self%number_of_modes
            idx(i) = i
        end do
        self%m = m_arr
        self%n = n_arr
        self%B_mn = B_mn
        call sorted_idx(self%n, idx)
        self%m = self%m(idx)
        self%n = self%n(idx)
        self%B_mn = self%B_mn(idx)
        self%n_min = self%n(1)
        call sorted_idx(self%m, idx)
        self%m = self%m(idx)
        self%n = self%n(idx)
        self%B_mn = self%B_mn(idx)

        ! print *, "m unsorted: ", m_arr
        ! print *, "n unsorted: ", n_arr
      print *, "Initialized Fourier field with ", self%number_of_modes, " sorted modes."
        ! print *, "m: ", self%m(1:min(self%number_of_modes,10)), "..."
        ! print *, "n: ", self%n(1:min(self%number_of_modes,10)), "..."

    end subroutine horner_fourier_field_init

    subroutine sorted_idx(arr_in, idx)
        integer, intent(in) :: arr_in(:)
        integer, intent(out) :: idx(:)

        integer, dimension(size(arr_in)) :: arr

        integer :: i, j, temp, temp_idx

        arr = arr_in

        do i = 1, size(arr)
            idx(i) = i
        end do

        do i = 1, size(arr) - 1
            do j = 1, size(arr) - i
                if (arr(j) > arr(j + 1)) then
                    temp = arr(j)
                    temp_idx = idx(j)
                    arr(j) = arr(j + 1)
                    idx(j) = idx(j + 1)
                    arr(j + 1) = temp
                    idx(j + 1) = temp_idx
                end if
            end do
        end do
    end subroutine sorted_idx

    subroutine compute_B_sqrtg_dB_dx(self, theta, phi, B_mod, sqrtg, dB_dx)
        class(horner_fourier_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod, sqrtg, dB_dx(3)

        print *, "horner_fourier_field_t does not provide sqrtg. "// &
            "Use compute_B_mod or compute_B_and_dB_dx instead!"
        error stop

    end subroutine compute_B_sqrtg_dB_dx

    subroutine compute_B_and_dB_dx(self, theta, phi, B_mod, dB_dx)
        class(horner_fourier_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod, dB_dx(3)

        integer :: k, idx_m, idx_n, prev_m, prev_n
        integer :: dm, dn, dn_min
        complex(dp) :: fac_m, fac_n, fac_n_min, i, prod, prod_m, prod_n
        complex(dp) :: B_temp, dB_dx_temp(3)

        fac_m = exp(cmplx(0.0_dp, theta, kind=dp))
        fac_n = exp(cmplx(0.0_dp, -phi*real(self%nfp, dp), kind=dp))
        fac_n_min = exp(cmplx(0.0_dp, -phi*real(self%n_min*self%nfp, dp), kind=dp))
        i = cmplx(0.0_dp, 1.0_dp, kind=dp)
        prod_m = cmplx(1.0_dp, 0.0_dp, kind=dp)
        prod_n = fac_n_min

        B_temp = 0.0_dp
        dB_dx_temp(1) = 0.0_dp
        dB_dx_temp(2) = 0.0_dp
        dB_dx_temp(3) = 0.0_dp

        k = 0
        prev_m = 0
        prev_n = self%n_min
        do k = 1, self%number_of_modes
            dm = self%m(k) - prev_m
            do idx_m = 1, dm
                prod_m = prod_m*fac_m
            end do
            dn = self%n(k) - prev_n
            prod_n = merge(prod_n, fac_n_min, dn >= 0)
            dn_min = self%n(k) - self%n_min
            dn = merge(dn, dn_min, dn >= 0)
            do idx_n = 1, dn
                prod_n = prod_n*fac_n
            end do
            prod = prod_m*prod_n
            B_temp = B_temp + self%B_mn(k)*prod
            dB_dx_temp(2) = dB_dx_temp(2) - self%B_mn(k)*real(self%m(k), dp)*prod
          dB_dx_temp(3) = dB_dx_temp(3) + self%B_mn(k)*real(self%n(k)*self%nfp, dp)*prod
            prev_m = self%m(k)
            prev_n = self%n(k)
        end do
        B_mod = real(B_temp, dp)
        dB_dx(1) = 0.0_dp
        dB_dx(2) = aimag(dB_dx_temp(2))
        dB_dx(3) = aimag(dB_dx_temp(3))
    end subroutine compute_B_and_dB_dx

    subroutine compute_B_mod(self, theta, phi, B_mod)
        class(horner_fourier_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod

        real(dp) :: dummy(3)

        call self%compute_B_and_dB_dx(theta, phi, B_mod, dummy)

    end subroutine compute_B_mod

    subroutine compute_nabla_s(self, theta, phi, nabla_s)
        class(horner_fourier_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: nabla_s

        nabla_s = 1.0_dp

    end subroutine compute_nabla_s

    real(dp) function rel_accuracy_B(self)
        class(horner_fourier_field_t), intent(in) :: self

        rel_accuracy_B = epsilon(1.0_dp)
    end function rel_accuracy_B

    subroutine get_covariant_components(self, B_theta_covariant, B_phi_covariant)
        class(horner_fourier_field_t), intent(in) :: self
        real(dp), intent(out) :: B_theta_covariant, B_phi_covariant

        print *, "horner_fourier_field_t does not provide covariant components."
        error stop

    end subroutine get_covariant_components

    !> Convert Fourier coefficients to a magnetic field on a grid using
    !! recurrence relations.
    subroutine ift_modes_to_B(m, n, B_mn, B, theta, phi, nfp, n_grid)
        use utils, only: linspace
        use constants, only: pi
        integer, intent(in) :: m(:)
            !! poloidal mode numbers (flat array, length mn_max)
        integer, intent(in) :: n(:)
            !! toroidal mode numbers normalised to nfp (flat array, length mn_max)
        real(dp), intent(in) :: B_mn(:)
            !! Fourier coefficients of B in Tesla (flat array, length mn_max)
        integer, intent(in), optional :: nfp
            !! number of field periods
        real(dp), intent(out), allocatable :: theta(:)
            !! poloidal angles in radians (flat array, length n_grid)
        real(dp), intent(out), allocatable :: phi(:)
            !! toroidal angles in radians (flat array, length n_grid)
        real(dp), intent(out), allocatable :: B(:, :)
            !! magnetic field in Tesla (2D array, size n_grid x n_grid)
        integer, intent(in), optional :: n_grid
            !! spline grid points per angle direction

        integer, parameter :: n_grid_default = 200

        integer :: n_theta, n_phi, i_theta, i_phi
        integer :: number_of_modes, i_mode
        complex(dp) :: imag
        complex(dp), dimension(:), allocatable :: exp_dtheta, exp_dphi
        complex(dp), dimension(:), allocatable :: exp_phi, exp_phase
        real(dp) :: dtheta, dphi, acc
        real(dp), allocatable :: m_real(:), n_real(:)
        real(dp) :: nfp_real

        if (.not. present(nfp)) then
            nfp_real = 1.0_dp
        elseif (nfp > 0) then
            nfp_real = real(nfp, dp)
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

        allocate (theta(n_theta), phi(n_phi))
        allocate (B(n_theta, n_phi))

        call linspace(0.0_dp, 2.0_dp*pi, n_theta, theta)
        call linspace(0.0_dp, 2.0_dp*pi/nfp_real, n_phi, phi)

        number_of_modes = size(m)
        allocate (m_real(number_of_modes), n_real(number_of_modes))
        m_real = real(m, dp)
        n_real = real(n, dp)*nfp_real

        allocate (exp_dtheta(number_of_modes), exp_dphi(number_of_modes))
        imag = (0.0_dp, 1.0_dp)
        dtheta = theta(2) - theta(1)
        exp_dtheta = exp(imag*m_real*dtheta)
        dphi = phi(2) - phi(1)
        exp_dphi = exp(-imag*n_real*dphi)

        allocate (exp_phi(number_of_modes), exp_phase(number_of_modes))
        exp_phi = (1.0_dp, 0.0_dp)

        do i_phi = 1, n_phi
            exp_phase = exp_phi
            do i_theta = 1, n_theta
                acc = 0.0_dp
                do i_mode = 1, number_of_modes
                    acc = acc + B_mn(i_mode)*real(exp_phase(i_mode), dp)
                end do
                B(i_theta, i_phi) = acc
                exp_phase = exp_phase*exp_dtheta
            end do
            exp_phi = exp_phi*exp_dphi
        end do
    end subroutine ift_modes_to_B

end module horner_fourier_field
