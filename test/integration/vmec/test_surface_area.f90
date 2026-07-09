!> Test that the surface area computed from the Boozer field matches
!> reference values obtained from simsopt (python/compute_surface_area.py).
program test_surface_area
    use constants, only: dp, pi
    use boozer_field, only: boozer_field_t
    use utils, only: not_same

    implicit none

    real(dp), parameter :: reltol = 1e-3, abstol = 1e-10
    character(len=*), parameter :: nc_filename = &
                      "input/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"

    type(boozer_field_t) :: bfield

    integer, parameter :: n_cases = 5
    real(dp) :: stor(n_cases)
    real(dp) :: area_ref(n_cases)

    real(dp) :: area
    integer :: case
    logical :: test_failed

    integer, parameter :: n_theta = 128, n_phi = 128
    real(dp) :: theta, phi, dtheta, dphi
    real(dp) :: nabla_s, sqrtg
    real(dp) :: dummy(4)
    integer :: i, j

    call bfield%init_from_vmec(nc_filename, &
                                  radial_spline_order=5, &
                                  angular_spline_order=5, &
                                  grid_refinement=3)

    test_failed = .false.
    stor = [0.1_dp, 0.3_dp, 0.5_dp, 0.7_dp, 0.9_dp]

    !> reference values from simsopt
    area_ref = [2.7980096435492493e+02_dp, &
                4.8598846172119789e+02_dp, &
                6.2887684715927901e+02_dp, &
                7.4580041608726401e+02_dp, &
                8.4759417849510817e+02_dp]

    dtheta = 2.0_dp*pi/real(n_theta, dp)
    dphi = 2.0_dp*pi/(bfield%nfp*real(n_phi, dp))

    do case = 1, n_cases
        call bfield%fix_to_surface(stor(case))
        area = 0.0_dp
        do j = 1, n_phi
            phi = (real(j, dp) - 0.5_dp)*dphi
            do i = 1, n_theta
                theta = (real(i, dp) - 0.5_dp)*dtheta
                call bfield%compute_nabla_s(theta, phi, nabla_s)
                call bfield%compute_B_sqrtg_dB_dx(theta, phi, &
                                                  dummy(1), &
                                                  sqrtg, &
                                                  dummy(2:4))
                area = area + abs(sqrtg*nabla_s)*dtheta*dphi
            end do
        end do
        area = area*bfield%nfp

        if (not_same(area, area_ref(case), &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "---------------------------------------------------"
            print *, "surface area mismatch at s =", stor(case)
            print *, "  computed: ", area
            print *, "  simsopt:  ", area_ref(case)
            print *, "  Relative error: ", abs(area - area_ref(case)) &
                /abs(area_ref(case))
            print *, "  Absolute error: ", abs(area - area_ref(case))
            test_failed = .true.
        end if
    end do

    if (test_failed) error stop

end program test_surface_area
