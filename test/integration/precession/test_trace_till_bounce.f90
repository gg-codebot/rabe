program test_trace_till_bounce
    use constants, only: dp, pi
    use boozer_field, only: boozer_field_t
    use utils, only: not_same
    use field_instance, only: initialize_field_instance
    use bounce, only: trace_orbit_till_bounce
    use params, only: params_init

    implicit none

    real(dp), parameter :: reltol = 1e-6, abstol = 1e-20
    real(dp), parameter :: abstol_for_zero = 1e-4
    character(len=*), parameter :: nc_filename = &
                      "input/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"

    type(boozer_field_t) :: bfield

    integer, parameter :: ndim_simple = 5
    integer, parameter :: ndim = 7
    real(dp), dimension(ndim) :: z_start, z_end
    real(dp) :: dtaumin_in
    integer, parameter :: n_bounces = 3
    integer :: i_bounce
    real(dp), dimension(n_bounces, ndim_simple) :: z_bounce_ref
    logical :: test_failed

    call bfield%init_from_vmec(nc_filename, &
                               radial_spline_order=5, &
                               angular_spline_order=5, &
                               grid_refinement=5, &
                               use_B_r_covariant=.false.)

    test_failed = .false.

    call initialize_field_instance(bfield)
    dtaumin_in = 2.0_dp*pi*bfield%R*1e2/128.0_dp
    call params_init(nfperiods=bfield%nfp, dtaumin_in=dtaumin_in)
    z_start = [0.5_dp, 6.2319269065015011_dp, 0.34860621153343652_dp, 1.0_dp, -0.1_dp, 0.0_dp, 0.0_dp]

    z_bounce_ref(1, :) = [4.0419130794376634e-01, 6.0365708273104213e+00, &
                          5.0430593251401215e+00, 1.0000000000000002e+00, &
                          -2.7822181191809768e-21]
    z_bounce_ref(2, :) = [4.0493842511538775e-01, 6.0434631646087142e+00, &
                          1.7569417065699922e+00, 9.9999999999999989e-01, &
                          2.1000115451565121e-21]
    z_bounce_ref(3, :) = [4.0164172878932580e-01, 6.0610633762317665e+00, &
                          4.8167018386562095e+00, 1.0000000000000000e+00, &
                          -5.2359841756454391e-21]

    do i_bounce = 1, n_bounces
        call trace_orbit_till_bounce(z_start, z_end)
        if (not_same(z_end(1:ndim_simple), &
                     z_bounce_ref(i_bounce, :), &
                     reltol, abstol)) then
            print *, "Test failed at bounce ", i_bounce
            print *, "Expected: ", z_bounce_ref(i_bounce, :)
            print *, "Got:      ", z_end(1:ndim_simple)
            print *, "Relative error z(1:4): ", abs((z_end(1:4) - &
                                                     z_bounce_ref(i_bounce, 1:4)) &
                                                    /z_bounce_ref(i_bounce, 1:4))
            print *, "Absolute error z(5): ", abs(z_end(5) - z_bounce_ref(i_bounce, 5))
            test_failed = .true.
        end if
        z_start = z_end
        z_start(ndim_simple + 1:ndim) = 0.0_dp
    end do

    if (test_failed) error stop
    print *, "--------------------------------------------------------"
    print *, "Test passed successfully."
    print *, "--------------------------------------------------------"

end program test_trace_till_bounce
