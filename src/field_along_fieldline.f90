module field_along_fieldline
    use constants, only: dp
    use field_base, only: field_t
    use fieldline_mod, only: fieldline_t

    implicit none
    private

    class(field_t), pointer, save :: field => null()
    type(fieldline_t), pointer, save :: fieldline => null()

    public :: B_mod_along_fieldline, dB_dphi_along_fieldline
    public :: set_field_and_fieldline, unset_field_and_fieldline

contains

    subroutine B_mod_along_fieldline(phi, B_mod)
        real(dp), dimension(:), intent(in) :: phi
        real(dp), dimension(:), intent(out) :: B_mod

        real(dp), dimension(size(phi)) :: theta
        integer :: idx

        theta = fieldline%get_theta(phi)
        do idx = 1, size(phi)
            call field%compute_B_mod(theta(idx), phi(idx), B_mod(idx))
        end do
    end subroutine B_mod_along_fieldline

    subroutine dB_dphi_along_fieldline(phi, dB_dphi)
        real(dp), dimension(:), intent(in) :: phi
        real(dp), dimension(:), intent(out) :: dB_dphi

        real(dp), dimension(size(phi)) :: theta
        real(dp) :: B_mod, dB_dx(3)
        integer :: idx

        theta = fieldline%get_theta(phi)
        do idx = 1, size(phi)
            call field%compute_B_and_dB_dx(theta(idx), phi(idx), &
                                           B_mod, dB_dx)
            dB_dphi(idx) = dB_dx(3) + fieldline%iota*dB_dx(2)
        end do
    end subroutine dB_dphi_along_fieldline

    subroutine set_field_and_fieldline(field_in, fieldline_in)
        class(field_t), intent(in), target :: field_in
        type(fieldline_t), intent(in), target :: fieldline_in

        call unset_field_and_fieldline()
        field => field_in
        fieldline => fieldline_in
    end subroutine set_field_and_fieldline

    subroutine unset_field_and_fieldline()
        field => null()
        fieldline => null()
    end subroutine unset_field_and_fieldline

end module field_along_fieldline
