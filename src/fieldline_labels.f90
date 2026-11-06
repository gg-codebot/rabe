module fieldline_labels
    use constants, only: dp, pi

    implicit none
    private

    type :: modes_t
        real(dp), dimension(:), allocatable :: cos_coeffs, sin_coeffs
        real(dp), dimension(:), allocatable :: mode_numbers
    end type modes_t

    type :: fieldline_modes_t
        type(modes_t) :: radial_drift
        type(modes_t) :: delta_eta
        type(modes_t) :: delta_aspect_ratio
    end type fieldline_modes_t

    public :: modes_t, fieldline_modes_t
    public :: fourier_transform_over_label, allocate_modes
    public :: get_labels, calc_iota_p

contains

    subroutine fourier_transform_over_label(flock, fieldline_modes)
        use fourier, only: real_ft
        use fieldline_mod, only: flock_of_fieldlines_t

        type(flock_of_fieldlines_t), intent(in) :: flock
        type(fieldline_modes_t), intent(out) :: fieldline_modes

        integer :: n_modes

        real(dp), dimension(size(flock%fieldlines)) :: shifted_label

        n_modes = size(flock%fieldlines)/2 + 1

        call allocate_modes(fieldline_modes%radial_drift, n_modes)
        call allocate_modes(fieldline_modes%delta_aspect_ratio, n_modes)
        call allocate_modes(fieldline_modes%delta_eta, n_modes)

        call real_ft(flock%fieldlines%xi_0, &
                     flock%fieldlines%radial_drift, &
                     fieldline_modes%radial_drift%cos_coeffs, &
                     fieldline_modes%radial_drift%sin_coeffs)

        call real_ft(flock%fieldlines%xi_0, &
                     flock%fieldlines%delta_aspect_ratio, &
                     fieldline_modes%delta_aspect_ratio%cos_coeffs, &
                     fieldline_modes%delta_aspect_ratio%sin_coeffs)

        shifted_label = flock%fieldlines%xi_0 - flock%iota_p
        call real_ft(shifted_label, &
                     flock%fieldlines%delta_eta, &
                     fieldline_modes%delta_eta%cos_coeffs, &
                     fieldline_modes%delta_eta%sin_coeffs)

    end subroutine fourier_transform_over_label

    subroutine allocate_modes(modes, n_modes)
        integer, intent(in) :: n_modes
        type(modes_t), intent(out) :: modes

        integer :: j

        allocate (modes%cos_coeffs(n_modes))
        allocate (modes%sin_coeffs(n_modes))
        allocate (modes%mode_numbers(n_modes))

        modes%cos_coeffs = 0.0_dp
        modes%sin_coeffs = 0.0_dp

        do j = 0, n_modes - 1
            modes%mode_numbers(j + 1) = real(j, kind=dp)
        end do
    end subroutine allocate_modes

    subroutine get_labels(max_n_fieldlines, iota, M_pol, N_tor, nfp, &
                          xi_0, approx_iota)
        use diophantine, only: rational_approx
        use utils, only: linspace
        integer, intent(in) :: max_n_fieldlines
        real(dp), intent(in) :: iota
        real(dp), intent(in) :: M_pol, N_tor, nfp
        real(dp), dimension(:), allocatable, intent(out) :: xi_0
        real(dp), intent(out) :: approx_iota

        real(dp) :: iota_p, iota_p_approx
        integer :: p, q
        integer :: n_fieldlines

        iota_p = calc_iota_p(iota, M_pol, N_tor, nfp)
        ! The symmetry points xi_0=pi and xi_0=pi-iota_p have to be either
        ! part of the label grid or lie symmetric between two labels. Only then
        ! are the sampled points of a symmetric function themselvese symmetric
        ! in respect to those points. The grid is automatic symmetric in respect
        ! to pi if its equidistant between 0 and 2pi. Additionally, if one finds an
        ! approximated iota so that iota_p/2pi is rational p/q, and takes q as
        ! number of points excluding the endpoint 2pi, then
        ! $$
        ! dxi_0 = 2pi/q = iota_p/p
        ! $$
        ! and pi-iota_p is a whole number of steps away from pi and therefore
        ! also either part of the grid or symmetric between two labels.
        call rational_approx(iota_p/(2.0_dp*pi), max_n_fieldlines, p, q)
        iota_p_approx = 2.0_dp*pi*p/q
        approx_iota = calc_iota(iota_p_approx, M_pol, N_tor, nfp)

        ! Any multiple of q would work, but we want to use the largest one that
        ! is smaller than max_n_fieldlines to get the best resolution.
        n_fieldlines = q*(max_n_fieldlines/q)
        allocate (xi_0(n_fieldlines))

        call linspace(0.0_dp, 2.0_dp*pi, n_fieldlines, &
                      xi_0, include_endpoint=.false.)

    end subroutine get_labels

    function calc_iota_p(iota, M_pol, N_tor, nfp) result(iota_p)
        real(dp), intent(in) :: iota
        real(dp), intent(in) :: M_pol, N_tor, nfp
        real(dp) :: iota_p

        if (abs(iota*M_pol - N_tor) < 1.0e-10_dp) then
            error stop "iota*M_pol - N_tor must not be zero (resonant)!"
        end if
        if (M_pol**2.0_dp + N_tor**2.0_dp < 1.0e-10_dp) then
            error stop "M_pol**2 + N_tor**2 must not be zero!"
        end if
        iota_p = sign(pi, iota*M_pol - N_tor)/(M_pol**2.0_dp + N_tor**2.0_dp)* &
                 (M_pol + &
                  nfp*(N_tor*iota + M_pol)/(iota*M_pol - N_tor))

    end function calc_iota_p

    function calc_iota(iota_p, M_pol, N_tor, nfp) result(iota)
        real(dp), intent(in) :: iota_p
        real(dp), intent(in) :: M_pol, N_tor, nfp
        real(dp) :: iota

        real(dp) :: z, sigma

        sigma = 1.0_dp
        z = iota_p*(M_pol**2.0_dp + N_tor**2.0_dp)/pi*sigma - M_pol
        iota = (M_pol*nfp + N_tor*z)/(M_pol*z - N_tor*nfp)

        if (iota*M_pol - N_tor < 0.0_dp) then
            sigma = -1.0_dp
            z = iota_p*(M_pol**2.0_dp + N_tor**2.0_dp)/pi*sigma - M_pol
            iota = (M_pol*nfp + N_tor*z)/(M_pol*z - N_tor*nfp)
        end if

    end function calc_iota

    function is_multiple_of_2pi(angle, tol_in)
        real(dp), intent(in) :: angle
        real(dp), intent(in), optional :: tol_in
        logical :: is_multiple_of_2pi

        real(dp) :: tol
        real(dp) :: remainder
        if (present(tol_in)) then
            tol = tol_in
        else
            tol = 3.0_dp*1e-2
        end if

        remainder = abs(mod(angle, 2.0_dp*pi))
        is_multiple_of_2pi = remainder < tol .or. abs(remainder - 2.0*pi) < tol
    end function is_multiple_of_2pi

end module fieldline_labels
