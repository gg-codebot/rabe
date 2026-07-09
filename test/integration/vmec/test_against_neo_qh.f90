!> Test that boozer_field_t produces the same results as the NEO code
!> for a quasi-helical equilibrium.
!> The reference values were obtained from a field evalution of a .bc file
!> converted from the original .nc file
program test_against_neo_qh
    use constants, only: dp
    use boozer_field, only: boozer_field_t
    use utils, only: not_same

    implicit none

    real(dp), parameter :: reltol = 1.2e-2, abstol = 1e-10
    real(dp), parameter :: abstol_for_zero = 1e-3
    ! Symmetry-zero dB_dx components (exact 0 in the NEO reference) come out
    ! as spline roundoff, measured <= 1.6e-10 with an optimized libneo build.
    real(dp), parameter :: abstol_roundoff = 1e-9
    character(len=*), parameter :: nc_filename = &
                      "input/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc"

    type(boozer_field_t) :: bfield

    integer, parameter :: n_cases = 5
    real(dp) :: stor(n_cases), theta(n_cases), phi(n_cases)

    real(dp) :: bmod_b, sqrtg_b, dB_dx_b(3), nabla_s
    real(dp) :: iota_b, B_theta_b, B_phi_b

    ! Neo reference values (from generate_neo_ref_qh)
    real(dp), parameter :: cm2m = 1e-2_dp
    real(dp), parameter :: psi_tor_edge_ref = 6.66284343900593610e+00_dp
    real(dp), parameter :: nfp_ref = 4.0_dp

    real(dp) :: iota_ref(n_cases)
    real(dp) :: B_theta_ref(n_cases), B_phi_ref(n_cases)
    real(dp) :: bmod_ref(n_cases), sqrtg_ref(n_cases)
    real(dp) :: dB_dx_ref(n_cases, 3)
    real(dp) :: nabla_s_ref(n_cases)

    integer :: case
    logical :: test_failed

    call bfield%init_from_vmec(nc_filename, &
                               radial_spline_order=5, &
                               angular_spline_order=5, &
                               grid_refinement=3)

    test_failed = .false.
    stor = [0.1_dp, 0.3_dp, 0.5_dp, 0.7_dp, 0.9_dp]
    theta = [0.0_dp, 1.0_dp, 3.14_dp, 0.5_dp, 2.0_dp]
    phi = [0.0_dp, 0.5_dp, 1.57_dp, 2.0_dp, 3.0_dp]

    iota_ref = [-1.24071689922512540e+00_dp, &
                -1.24292805961154729e+00_dp, &
                -1.24354459000000217e+00_dp, &
                -1.24415222615299226e+00_dp, &
                -1.24472774238178752e+00_dp]

    B_theta_ref = [1.95255977851786137e-16_dp, &
                   5.86292403810843430e-16_dp, &
                   -6.21031004000000052e-16_dp, &
                   -4.40056927220479126e-16_dp, &
                   -2.19425745400124919e-17_dp]

    B_phi_ref = [1.02238518141380709e+02_dp, &
                 1.02238488245105472e+02_dp, &
                 1.02238467200000187e+02_dp, &
                 1.02238454689203053e+02_dp, &
                 1.02238446704235940e+02_dp]

    bmod_ref = [5.58389078029500130e+00_dp, &
                6.44701243605253538e+00_dp, &
                6.64030869700747317e+00_dp, &
                6.43027846079623178e+00_dp, &
                5.80708194880627726e+00_dp]

    sqrtg_ref = [-2.18474283967293203e+07_dp, &
                 -1.63891719501446038e+07_dp, &
                 -1.54488942118170653e+07_dp, &
                 -1.64745789734431263e+07_dp, &
                 -2.02003056883056164e+07_dp]
    sqrtg_ref = sqrtg_ref*(cm2m**3.0_dp)

    dB_dx_ref(1, :) = [-1.44077730567127604e+00_dp, &
                       0.00000000000000000e+00_dp, &
                       0.00000000000000000e+00_dp]
    dB_dx_ref(2, :) = [1.01155913223093696e+00_dp, &
                       8.13915228819103859e-02_dp, &
                       3.29782295582588203e-01_dp]
    dB_dx_ref(3, :) = [8.24117782542400312e-01_dp, &
                       3.85688857474143776e-03_dp, &
                       1.52445741914789129e-02_dp]
    dB_dx_ref(4, :) = [4.34149361141258405e-01_dp, &
                       6.94592020903635365e-01_dp, &
                       2.78307702440704308e+00_dp]
    dB_dx_ref(5, :) = [-2.83188161564565527e-02_dp, &
                       9.01242770896047984e-01_dp, &
                       3.60496554369149980e+00_dp]

    nabla_s_ref = [0.75369886165216027_dp, &
                   0.67255511319280181_dp, &
                   2.0144853453318401_dp, &
                   0.80518369474188856_dp, &
                   0.77957223538352005_dp]

    ! Compare global quantities
    if (not_same(bfield%psi_tor_edge, psi_tor_edge_ref, &
                 reltol_in=reltol, abstol_in=abstol)) then
        print *, "psi_tor_edge mismatch:"
        print *, "  boozer: ", bfield%psi_tor_edge
        print *, "  neo:    ", psi_tor_edge_ref
        print *, "Relative error: ", abs(bfield%psi_tor_edge - psi_tor_edge_ref) &
            /abs(psi_tor_edge_ref)
        print *, "Absolute error: ", abs(bfield%psi_tor_edge - psi_tor_edge_ref)
        test_failed = .true.
    end if

    if (not_same(bfield%nfp, nfp_ref, &
                 reltol_in=reltol, abstol_in=abstol)) then
        print *, "nfp mismatch:"
        print *, "  boozer: ", bfield%nfp
        print *, "  neo:    ", nfp_ref
        print *, "Relative error: ", abs(bfield%nfp - nfp_ref)/abs(nfp_ref)
        print *, "Absolute error: ", abs(bfield%nfp - nfp_ref)
        test_failed = .true.
    end if

    do case = 1, n_cases
        call bfield%fix_to_surface(stor(case))

        call bfield%get_iota(stor(case), iota_b)
        call bfield%get_covariant_components(B_theta_b, B_phi_b)

        if (not_same(iota_b, iota_ref(case), &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "iota mismatch at case ", case
            print *, "  boozer: ", iota_b
            print *, "  neo:    ", iota_ref(case)
            print *, "Relative error: ", abs(iota_b - iota_ref(case)) &
                /abs(iota_ref(case))
            print *, "Absolute error: ", abs(iota_b - iota_ref(case))
            test_failed = .true.
        end if

        if (not_same(B_theta_b, B_theta_ref(case), &
                     reltol_in=0.0_dp, abstol_in=abstol_for_zero)) then
            print *, "B_theta_covariant mismatch at case ", case
            print *, "  boozer: ", B_theta_b
            print *, "  neo:    ", B_theta_ref(case)
            print *, "Relative error: ", abs(B_theta_b - B_theta_ref(case)) &
                /abs(B_theta_ref(case))
            print *, "Absolute error: ", abs(B_theta_b - B_theta_ref(case))
            test_failed = .true.
        end if

        if (not_same(B_phi_b, B_phi_ref(case), &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "B_phi_covariant mismatch at case ", case
            print *, "  boozer: ", B_phi_b
            print *, "  neo:    ", B_phi_ref(case)
            print *, "Relative error: ", abs(B_phi_b - B_phi_ref(case)) &
                /abs(B_phi_ref(case))
            print *, "Absolute error: ", abs(B_phi_b - B_phi_ref(case))
            test_failed = .true.
        end if

        call bfield%compute_B_sqrtg_dB_dx(theta(case), phi(case), &
                                          bmod_b, sqrtg_b, dB_dx_b)

        if (not_same(bmod_b, bmod_ref(case), &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "B_mod mismatch at case ", case
            print *, "  boozer: ", bmod_b
            print *, "  neo:    ", bmod_ref(case)
            print *, "Relative error: ", abs(bmod_b - bmod_ref(case)) &
                /abs(bmod_ref(case))
            print *, "Absolute error: ", abs(bmod_b - bmod_ref(case))
            test_failed = .true.
        end if

        if (not_same(sqrtg_b, sqrtg_ref(case), &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "sqrtg mismatch at case ", case
            print *, "  boozer: ", sqrtg_b
            print *, "  neo:    ", sqrtg_ref(case)
            print *, "Relative error: ", abs(sqrtg_b - sqrtg_ref(case)) &
                /abs(sqrtg_ref(case))
            print *, "Absolute error: ", abs(sqrtg_b - sqrtg_ref(case))
            test_failed = .true.
        end if

        if (not_same(dB_dx_b, dB_dx_ref(case, :), &
                     reltol_in=reltol, abstol_in=abstol_roundoff)) then
            print *, "dB_dx mismatch at case ", case
            print *, "  boozer: ", dB_dx_b
            print *, "  neo:    ", dB_dx_ref(case, :)
            print *, "Relative error: ", abs(dB_dx_b - dB_dx_ref(case, :)) &
                /abs(dB_dx_ref(case, :))
            print *, "Absolute error: ", abs(dB_dx_b - dB_dx_ref(case, :))
            test_failed = .true.
        end if

        call bfield%compute_nabla_s(theta(case), phi(case), nabla_s)

        if (not_same(nabla_s, nabla_s_ref(case), &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "nabla_s mismatch at case ", case
            print *, "  boozer: ", nabla_s
            print *, "  neo:    ", nabla_s_ref(case)
            print *, "Relative error: ", abs(nabla_s - nabla_s_ref(case)) &
                /abs(nabla_s_ref(case))
            print *, "Absolute error: ", abs(nabla_s - nabla_s_ref(case))
            test_failed = .true.
        end if
    end do

    if (test_failed) error stop

end program test_against_neo_qh
