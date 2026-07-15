module mock_field_3d
    use constants, only: dp, pi
    use field_base, only: field_t
    use field_base, only: field_3D_t

    implicit none

    type, extends(field_3D_t) :: mock_field_3d_t
        class(field_t), allocatable :: field_2D
        real(dp) :: psi_tor_edge = 1.0_dp
        real(dp) :: iota = 1.0_dp
    contains
        procedure :: mock_field_3d_init
        procedure :: compute_B_sqrtg_dB_dx
        procedure :: compute_B_and_dB_dx
        procedure :: compute_B_mod
        procedure :: compute_nabla_s
        procedure :: rel_accuracy_B
        procedure :: get_covariant_components
        procedure :: evaluate
        procedure :: get_iota
    end type mock_field_3d_t

contains

    subroutine mock_field_3d_init(self, field_2D_in)
        class(mock_field_3d_t), intent(out) :: self
        class(field_t), intent(in) :: field_2D_in

        allocate (self%field_2D, source=field_2D_in)

    end subroutine mock_field_3d_init

    subroutine compute_B_sqrtg_dB_dx(self, theta, phi, B_mod, sqrtg, dB_dx)
        class(mock_field_3d_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod, sqrtg, dB_dx(3)

        call self%field_2D%compute_B_sqrtg_dB_dx(theta, phi, B_mod, sqrtg, dB_dx)

    end subroutine compute_B_sqrtg_dB_dx

    subroutine compute_B_and_dB_dx(self, theta, phi, B_mod, dB_dx)
        class(mock_field_3d_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod, dB_dx(3)

        call self%field_2D%compute_B_and_dB_dx(theta, phi, B_mod, dB_dx)

    end subroutine compute_B_and_dB_dx

    subroutine compute_B_mod(self, theta, phi, B_mod)
        class(mock_field_3d_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod

        call self%field_2D%compute_B_mod(theta, phi, B_mod)

    end subroutine compute_B_mod

    subroutine compute_nabla_s(self, theta, phi, nabla_s)
        class(mock_field_3d_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: nabla_s

        call self%field_2D%compute_nabla_s(theta, phi, nabla_s)

    end subroutine compute_nabla_s

    real(dp) function rel_accuracy_B(self)
        class(mock_field_3d_t), intent(in) :: self

        rel_accuracy_B = self%field_2D%rel_accuracy_B()
    end function rel_accuracy_B

    subroutine get_covariant_components(self, B_theta_covariant, B_phi_covariant)
        class(mock_field_3d_t), intent(in) :: self
        real(dp), intent(out) :: B_theta_covariant, B_phi_covariant

        call self%field_2D%get_covariant_components(B_theta_covariant, &
                                                    B_phi_covariant)
    end subroutine get_covariant_components

    subroutine evaluate(self, x, bmod, sqrtg, bder, &
                        hcovar, hctrvr, hcurl)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
        use constants, only: machine_eps
        class(mock_field_3d_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: bmod, sqrtg
        real(dp), intent(out) :: bder(3), hcovar(3), hctrvr(3), hcurl(3)

        real(dp) :: B_vartheta, B_varphi
        real(dp) :: B_mod, dB_dx(3)

        if (any(ieee_is_nan(x))) then
            print *, "Input x is NaN!"
            print *, "x = ", x
            error stop
        end if

        call self%field_2D%compute_B_and_dB_dx(x(2), x(3), B_mod, dB_dx)

        bmod = B_mod
        bder = dB_dx/B_mod
        call self%get_covariant_components(B_vartheta, B_varphi)
        ! covariant unit vector consistent with covariant components of B
        hcovar(1) = 0.0_dp ! choice for mock field to have zero radial covariant
        hcovar(2) = B_vartheta/B_mod
        hcovar(3) = B_varphi/B_mod
        ! contravariant unit vector consistent with Boozer relation
        ! B_ctrvr_varphi = B_squared/(B_varphi + iota*B_vartheta)
        ! and defintion of iota
        hctrvr(1) = 0.0_dp
        hctrvr(3) = B_mod/(B_varphi + self%iota*B_vartheta)
        hctrvr(2) = self%iota*hctrvr(3)

        sqrtg = self%psi_tor_edge/(B_mod*hctrvr(3))

        ! as for the mock field dB/ds=0 and radial covariant component is zero
        hcurl(1) = -(hcovar(3)*dB_dx(2) - hcovar(2)*dB_dx(3))/B_mod
        hcurl(2) = 0.0_dp
        hcurl(3) = 0.0_dp
        hcurl = hcurl/sqrtg

        if (ieee_is_nan(bmod)) then
            print *, "bmod is NaN!"
            error stop
        end if
        if (ieee_is_nan(sqrtg)) then
            print *, "sqrtg is NaN!"
            error stop
        end if
        if (any(ieee_is_nan(bder))) then
            print *, "bder is NaN!"
            error stop
        end if
        if (any(ieee_is_nan(hcovar))) then
            print *, "hcovar is NaN!"
            error stop
        end if
        if (any(ieee_is_nan(hctrvr))) then
            print *, "hctrvr is NaN!"
            error stop
        end if
        if (any(ieee_is_nan(hcurl))) then
            print *, "hcurl is NaN!"
            error stop
        end if
        if (abs(bmod) < machine_eps) then
            print *, "bmod is zero!"
            error stop
        end if
        if (abs(sqrtg) < machine_eps) then
            print *, "sqrtg is zero!"
            error stop
        end if
    end subroutine evaluate

    subroutine get_iota(self, s_tor, iota)
        class(mock_field_3d_t), intent(in) :: self
        real(dp), intent(in) :: s_tor
        real(dp), intent(out) :: iota

        real(dp) :: dummy(11), hctrvr(3)
        real(dp) :: x(3)

        x(1) = s_tor
        x(2) = 0.0_dp
        x(3) = 0.0_dp

        call self%evaluate(x, dummy(1), dummy(2), dummy(3:5), &
                           dummy(6:8), hctrvr, dummy(9:11))

        iota = hctrvr(2)/hctrvr(3)

    end subroutine get_iota

end module mock_field_3d
