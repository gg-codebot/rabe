module fieldline_mod
    use constants, only: dp, pi

    implicit none
    private

    public :: fieldline_t, flock_of_fieldlines_t

    type :: flock_of_fieldlines_t
        !! Collection of field lines on a flux surface and surface properties.
        type(fieldline_t), allocatable :: fieldlines(:)
        real(dp) :: iota
        real(dp) :: M_pol
        real(dp) :: N_tor
        real(dp) :: nfp
        real(dp) :: iota_p
        real(dp) :: I_ref
        real(dp) :: eta_b
    end type flock_of_fieldlines_t

    type :: fieldline_t
        !! Single magnetic field line: starting label, extrema, and bounce integrals.
        real(dp) :: xi_0
        real(dp) :: theta_0
        real(dp) :: phi_0
        real(dp) :: iota
        real(dp) :: phi_max(2)
        real(dp) :: phi_max_error(2)
        real(dp) :: B_max(2)

        real(dp) :: eta_b
        real(dp) :: delta_eta
        real(dp) :: integral_lambda_b_over_B_squared
        real(dp) :: delta_aspect_ratio
        real(dp) :: integral_one_over_B_squared
        real(dp) :: radial_drift
        real(dp) :: integral_nabla_s_over_B_squared
    contains
        generic :: get_theta => get_theta_scalar, get_theta_array
        procedure, private :: get_theta_scalar
        procedure, private :: get_theta_array
    end type fieldline_t

contains

    function get_theta_scalar(self, phi) result(theta)
        class(fieldline_t), intent(in) :: self
        real(dp) :: phi

        real(dp) :: theta

        theta = (phi - self%phi_0)*self%iota + self%theta_0
    end function get_theta_scalar

    function get_theta_array(self, phi) result(theta)
        class(fieldline_t), intent(in) :: self
        real(dp), dimension(:) :: phi

        real(dp), dimension(size(phi)) :: theta

        theta = (phi - self%phi_0)*self%iota + self%theta_0
    end function get_theta_array

end module fieldline_mod
