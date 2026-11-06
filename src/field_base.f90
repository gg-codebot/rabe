module field_base
    use constants, only: dp
    implicit none
    private

    public :: field_t, field_3D_t

    type, abstract :: field_t
        !! Abstract base type for magnetic fields in Boozer coordinates.
        !!
        !! Concrete implementations supply B_mod, sqrt(g), gradients, and covariant components.
    contains
        procedure(compute_B_sqrtg_dB_dx), deferred :: compute_B_sqrtg_dB_dx
        procedure(compute_B_and_dB_dx), deferred :: compute_B_and_dB_dx
        procedure(compute_B_mod), deferred :: compute_B_mod
        procedure(compute_nabla_s), deferred :: compute_nabla_s
        procedure(rel_accuracy_B), deferred :: rel_accuracy_B
        procedure(get_covariant_components), deferred :: get_covariant_components
    end type field_t

    interface
        subroutine compute_B_sqrtg_dB_dx(self, theta, phi, B_mod, sqrtg, dB_dx)
            import :: field_t, dp
            class(field_t), intent(in) :: self
            real(dp), intent(in) :: theta, phi
            real(dp), intent(out) :: B_mod, sqrtg, dB_dx(3)
        end subroutine
    end interface

    interface
        subroutine compute_B_and_dB_dx(self, theta, phi, B_mod, dB_dx)
            import :: field_t, dp
            class(field_t), intent(in) :: self
            real(dp), intent(in) :: theta, phi
            real(dp), intent(out) :: B_mod, dB_dx(3)
        end subroutine
    end interface

    interface
        subroutine compute_B_mod(self, theta, phi, B_mod)
            import :: field_t, dp
            class(field_t), intent(in) :: self
            real(dp), intent(in) :: theta, phi
            real(dp), intent(out) :: B_mod
        end subroutine
    end interface

    interface
        subroutine compute_nabla_s(self, theta, phi, nabla_s)
            import :: field_t, dp
            class(field_t), intent(in) :: self
            real(dp), intent(in) :: theta, phi
            real(dp), intent(out) :: nabla_s
        end subroutine
    end interface

    interface
        !> Relative accuracy of B_mod evaluations for this field representation.
        real(dp) function rel_accuracy_B(self)
            import :: field_t, dp
            class(field_t), intent(in) :: self
        end function
    end interface

    interface
        !> Covariant components of B (flux-surface constants in Boozer coordinates).
        subroutine get_covariant_components(self, B_theta_covariant, B_phi_covariant)
            import :: field_t, dp
            class(field_t), intent(in) :: self
            real(dp), intent(out) :: B_theta_covariant, B_phi_covariant
        end subroutine
    end interface

    type, abstract, extends(field_t) :: field_3D_t
    contains
        procedure(evaluate), deferred :: evaluate
    end type field_3D_t

    interface
        subroutine evaluate(self, x, bmod, sqrtg, bder, &
                            hcovar, hctrvr, hcurl)
            import :: field_3D_t, dp
            class(field_3D_t), intent(in) :: self
            real(dp), intent(in) :: x(3)
            real(dp), intent(out) :: bmod, sqrtg
            real(dp), intent(out) :: bder(3), hcovar(3), hctrvr(3), hcurl(3)
        end subroutine evaluate
    end interface

end module field_base
