module shaped_well_field
    use constants, only: dp
    use field_base, only: field_t

    implicit none

    !> B = B_0*(1 + eps*cos(u)), u = chi + delta*cos(theta)*sin(chi),
    !> chi = M_pol*theta - N_tor*phi.
    !> Along every field line u is monotone (for small enough delta), so all lines
    !> share the same B_min and B_max, but the well width depends on the label.
    type, extends(field_t) :: shaped_well_field_t
        real(dp) :: M_pol, N_tor
        real(dp) :: B_0, eps, delta
    contains
        procedure :: shaped_well_field_init
        procedure :: compute_B_sqrtg_dB_dx
        procedure :: compute_B_and_dB_dx
        procedure :: compute_B_mod
        procedure :: compute_nabla_s
        procedure :: rel_accuracy_B
        procedure :: get_covariant_components
    end type shaped_well_field_t

contains

    subroutine shaped_well_field_init(self, M_pol, N_tor, B_0, eps, delta)
        class(shaped_well_field_t), intent(out) :: self
        real(dp), intent(in) :: M_pol, N_tor
        real(dp), intent(in) :: B_0, eps, delta

        self%M_pol = M_pol
        self%N_tor = N_tor
        self%B_0 = B_0
        self%eps = eps
        self%delta = delta
    end subroutine shaped_well_field_init

    subroutine compute_B_sqrtg_dB_dx(self, theta, phi, B_mod, sqrtg, dB_dx)
        class(shaped_well_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod, sqrtg, dB_dx(3)

        print *, "The type shaped_well_field does not provide sqrtg."// &
            "Use compute_B_mod or compute_B_and_dB_dx instead!"
        error stop

    end subroutine compute_B_sqrtg_dB_dx

    subroutine compute_B_and_dB_dx(self, theta, phi, B_mod, dB_dx)
        class(shaped_well_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod, dB_dx(3)

        real(dp) :: chi, u, du_dtheta, du_dphi

        chi = self%M_pol*theta - self%N_tor*phi
        u = chi + self%delta*cos(theta)*sin(chi)
        du_dtheta = self%M_pol*(1.0_dp + self%delta*cos(theta)*cos(chi)) &
                    - self%delta*sin(theta)*sin(chi)
        du_dphi = -self%N_tor*(1.0_dp + self%delta*cos(theta)*cos(chi))

        B_mod = self%B_0*(1.0_dp + self%eps*cos(u))
        dB_dx(1) = 0.0_dp
        dB_dx(2) = -self%B_0*self%eps*sin(u)*du_dtheta
        dB_dx(3) = -self%B_0*self%eps*sin(u)*du_dphi
    end subroutine compute_B_and_dB_dx

    subroutine compute_B_mod(self, theta, phi, B_mod)
        class(shaped_well_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod

        real(dp) :: chi, u

        chi = self%M_pol*theta - self%N_tor*phi
        u = chi + self%delta*cos(theta)*sin(chi)
        B_mod = self%B_0*(1.0_dp + self%eps*cos(u))
    end subroutine compute_B_mod

    subroutine compute_nabla_s(self, theta, phi, nabla_s)
        class(shaped_well_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: nabla_s

        nabla_s = 1.0_dp

    end subroutine compute_nabla_s

    real(dp) function rel_accuracy_B(self)
        class(shaped_well_field_t), intent(in) :: self

        rel_accuracy_B = 1e-14_dp
    end function rel_accuracy_B

    !> rabe works internally with modulus of B only
    !> choose covariants so that B^\varphi = B^2 (no-op for tests).
    subroutine get_covariant_components(self, B_theta_covariant, B_phi_covariant)
        class(shaped_well_field_t), intent(in) :: self
        real(dp), intent(out) :: B_theta_covariant, B_phi_covariant

        B_theta_covariant = 0.0_dp
        B_phi_covariant = 1.0_dp
    end subroutine get_covariant_components

end module shaped_well_field
