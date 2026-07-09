!> Test that boozer_field_t produces the same results as the NEO code
!> The reference values were obtained from a field evalution of a .bc file
!> converted from the original .nc file
program test_against_neo_qa
    use constants, only: dp
    use boozer_field, only: boozer_field_t
    use utils, only: not_same

    implicit none

    real(dp), parameter :: reltol = 1.4e-2, abstol = 1e-10
    real(dp), parameter :: abstol_for_zero = 1e-4
    character(len=*), parameter :: nc_filename = &
                      "input/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"

    type(boozer_field_t) :: bfield

    integer, parameter :: n_cases = 5
    real(dp) :: stor(n_cases), theta(n_cases), phi(n_cases)

    real(dp) :: bmod_b, sqrtg_b, dB_dx_b(3), nabla_s_b
    real(dp) :: iota_b, B_theta_b, B_phi_b

    ! Neo reference values (from neo_reference_values.dat)
    real(dp), parameter :: cm2m = 1e-2_dp
    real(dp), parameter :: psi_tor_edge_ref = 8.03450599209855731e+00_dp
    real(dp), parameter :: nfp_ref = 2.0_dp

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

    iota_ref = [4.22422154730570931e-01_dp, &
                4.20762476348729708e-01_dp, &
                4.19208746000000243e-01_dp, &
                4.17769165311363710e-01_dp, &
                4.16425934774603823e-01_dp]

    B_theta_ref = [3.02654553841740003e-16_dp, &
                   -5.30847580644957470e-17_dp, &
                   5.67254575999999982e-16_dp, &
                   3.88992237713905187e-18_dp, &
                   4.09913822581704287e-16_dp]

    B_phi_ref = [6.35881257119542482e+01_dp, &
                 6.35881130037971829e+01_dp, &
                 6.35881071999999961e+01_dp, &
                 6.35881041064855523e+01_dp, &
                 6.35881021533644812e+01_dp]

    bmod_ref = [5.70341902299086545e+00_dp, &
                5.70108512074173568e+00_dp, &
                6.37422408379196259e+00_dp, &
                5.42184282133529383e+00_dp, &
                6.10812813169421798e+00_dp]

    sqrtg_ref = [-1.57059656572236344e+07_dp, &
                 -1.57188245204609483e+07_dp, &
                 -1.25742017125101890e+07_dp, &
                 -1.73796573177368827e+07_dp, &
                 -1.36936355354960281e+07_dp]
    sqrtg_ref = sqrtg_ref*(cm2m**3.0_dp)

    dB_dx_ref(1, :) = [-1.02319025931563679e+00_dp, &
                       0.00000000000000000e+00_dp, &
                       0.00000000000000000e+00_dp]
    dB_dx_ref(2, :) = [-3.73649657269343738e-01_dp, &
                       2.77175396200966873e-01_dp, &
                       -1.30692773943976294e-03_dp]
    dB_dx_ref(3, :) = [4.69504795454586898e-01_dp, &
                       8.61561331571138438e-04_dp, &
                       -1.50754488577436443e-04_dp]
    dB_dx_ref(4, :) = [-3.54433115319706837e-01_dp, &
                       2.14487046410498761e-01_dp, &
                       -3.93159210289616753e-03_dp]
    dB_dx_ref(5, :) = [7.75627349470958266e-02_dp, &
                       6.17949257583749190e-01_dp, &
                       -5.95561440956370904e-03_dp]

    nabla_s_ref = [0.69190314914019782e+01_dp, &
                   0.78340849662291232e+01_dp, &
                   0.12852437144249926e+01_dp, &
                   0.53259697504179320e+01_dp, &
                   0.10693676025791250e+02_dp]/psi_tor_edge_ref
    nabla_s_ref = [0.75803529773215883_dp, &
                   1.0210847193060111_dp, &
                   0.61913969252584766_dp, &
                   0.72623895266343930_dp, &
                   1.0088316873581507_dp]

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

        ! Skip dB/dphi (component 3) — for QA it is at the spline noise floor
        if (not_same(dB_dx_b(1:2), dB_dx_ref(case, 1:2), &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "dB_dx mismatch at case ", case
            print *, "  boozer: ", dB_dx_b(1:2)
            print *, "  neo:    ", dB_dx_ref(case, 1:2)
            print *, "Relative error: ", abs(dB_dx_b(1:2) - dB_dx_ref(case, 1:2)) &
                /abs(dB_dx_ref(case, 1:2))
            print *, "Absolute error: ", abs(dB_dx_b(1:2) - dB_dx_ref(case, 1:2))
            test_failed = .true.
        end if

        call bfield%compute_nabla_s(theta(case), phi(case), nabla_s_b)

        if (not_same(nabla_s_b, nabla_s_ref(case), &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, "nabla_s mismatch at case ", case
            print *, "  boozer: ", nabla_s_b
            print *, "  neo:    ", nabla_s_ref(case)
            print *, "Relative error: ", abs(nabla_s_b - nabla_s_ref(case)) &
                /abs(nabla_s_ref(case))
            print *, "Absolute error: ", abs(nabla_s_b - nabla_s_ref(case))
            test_failed = .true.
        end if
    end do

    if (test_failed) error stop

end program test_against_neo_qa
