module surface_average_mod
    use constants, only: dp, pi, machine_eps
    implicit none
    private

    type :: surface_average_t
        !! flux surface averaged quantities computed by calc_surface_averages.
        real(dp) :: normalization
        real(dp) :: B_squared
        real(dp) :: lambda_b
        real(dp) :: nabla_s
    end type surface_average_t

    public :: surface_average_t, calc_surface_averages

contains

    !> Compute bounce-averaged surface quantities from a flock of field lines.
    !!
    !! Requires at least two field lines in flock.
    subroutine calc_surface_averages(flock, surface_average)
        use fieldline_mod, only: flock_of_fieldlines_t
        use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
        type(flock_of_fieldlines_t), intent(in) :: flock
        type(surface_average_t), intent(out) :: surface_average

        real(dp), dimension(size(flock%fieldlines)) :: well_lengths

        integer :: n_fieldlines
        real(dp) :: dxi_0

        n_fieldlines = size(flock%fieldlines)
        if (n_fieldlines < 2) then
   print *, "error: at least two fieldlines are required to calculate surface averages."
            error stop
        end if

        dxi_0 = (flock%fieldlines(n_fieldlines)%xi_0 - flock%fieldlines(1)%xi_0)/ &
                real(n_fieldlines - 1, kind=dp)

        well_lengths = flock%fieldlines%phi_max(2) - flock%fieldlines%phi_max(1)
        surface_average%normalization = sum( &
                                        flock%fieldlines%integral_one_over_B_squared &
                                        )*dxi_0
        if (abs(surface_average%normalization) < machine_eps) then
            print *, "error: surface average normalization is zero."
            error stop
        end if
        surface_average%B_squared = sum(well_lengths)*dxi_0/ &
                                    surface_average%normalization
        surface_average%lambda_b = sum( &
                                   flock%fieldlines%integral_lambda_b_over_B_squared &
                                   )*dxi_0/surface_average%normalization
        surface_average%nabla_s = sum( &
                                  flock%fieldlines%integral_nabla_s_over_B_squared &
                                  )*dxi_0/surface_average%normalization

        if (surface_average%lambda_b < machine_eps) then
            print *, "error: average lambda_b <= 0."
            error stop
        end if
        if (surface_average%B_squared < machine_eps) then
            print *, "error: average B_squared <= 0."
            error stop
        end if
        if (ieee_is_nan(surface_average%nabla_s)) then
            print *, "error: average sqrt_g11 is NaN."
            error stop
        end if
    end subroutine calc_surface_averages

end module surface_average_mod
