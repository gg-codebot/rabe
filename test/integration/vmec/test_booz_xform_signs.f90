!> Sign-convention test for the booz_xform loader path.
!>
!> Builds a synthetic boozmn dataset with positive and negative ixn modes,
!> writes it via libneo write_boozmn, loads it through init_from_boozmn, and
!> checks Bmod against the analytic Fourier sum with the booz_xform standard
!> cos(m*theta - n*zeta) phase.
program test_booz_xform_signs
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use constants, only: pi
    use boozmn_file, only: boozmn_data_t, write_boozmn
    use boozer_field, only: boozer_field_t
    use utils, only: not_same

    implicit none

    character(len=*), parameter :: boozmn_file_name = 'test_booz_xform_signs.nc'
    integer, parameter :: nmode = 4
    integer, parameter :: ns_full = 5
    integer, parameter :: nsurf = 3
    integer, dimension(nmode), parameter :: ixm = [0, 1, 1, 2]
    integer, dimension(nmode), parameter :: ixn = [0, 1, -1, 0]
    real(dp), dimension(nmode), parameter :: bmn = [2.0_dp, 0.1_dp, -0.07_dp, 0.02_dp]
    real(dp), parameter :: twopi = 2.0_dp*pi
    real(dp), parameter :: stor = 0.5_dp
    real(dp), parameter :: theta = twopi*7.0_dp/48.0_dp
    real(dp), parameter :: zeta = twopi*11.0_dp/96.0_dp
    real(dp), parameter :: reltol = 1.0e-6_dp
    real(dp), parameter :: abstol = 1.0e-10_dp

    type(boozmn_data_t) :: d
    type(boozer_field_t) :: field
    real(dp) :: bmod, expected

    d%ns = ns_full
    d%nfp = 1
    d%nmodes = nmode
    d%nsurf = nsurf
    d%asym = .false.
    allocate (d%jlist(nsurf), d%ixm(nmode), d%ixn(nmode))
    allocate (d%iota(ns_full), d%buco(ns_full), d%bvco(ns_full), d%phi(ns_full))
    allocate (d%bmnc(nmode, nsurf), d%rmnc(nmode, nsurf), &
              d%zmns(nmode, nsurf), d%pmns(nmode, nsurf))

    d%jlist = [2, 3, 4]
    d%ixm = ixm
    d%ixn = ixn
    d%iota = [0.5_dp, 0.55_dp, 0.6_dp, 0.65_dp, 0.7_dp]
    d%buco = [0.01_dp, 0.02_dp, 0.03_dp, 0.04_dp, 0.05_dp]
    d%bvco = [0.2_dp, 0.21_dp, 0.22_dp, 0.23_dp, 0.24_dp]
    d%phi = [0.0_dp, 0.2_dp, 0.5_dp, 0.9_dp, 1.3_dp]
    d%bmnc = spread(bmn, dim=2, ncopies=nsurf)
    d%rmnc = 0.0_dp
    d%rmnc(1, :) = [1.8_dp, 1.9_dp, 2.0_dp]
    d%zmns = 0.0_dp
    d%pmns = 0.0_dp

    call write_boozmn(boozmn_file_name, d)
    call field%init_from_boozmn(boozmn_file_name)
    if (.not. field%initialized) error stop 'booz_xform field not initialized'

    call field%fix_to_surface(stor)
    call field%compute_B_mod(theta, zeta, bmod)
    expected = expected_bmod(theta, zeta)

    if (not_same(bmod, expected, reltol_in=reltol, abstol_in=abstol)) then
        print *, 'booz_xform signed-mode Bmod:', bmod, ' expected:', expected
        error stop 'booz_xform signed-mode Bmod mismatch'
    end if

contains

    function expected_bmod(theta_in, zeta_in) result(value)
        real(dp), intent(in) :: theta_in, zeta_in
        real(dp) :: value
        integer :: imode

        value = 0.0_dp
        do imode = 1, nmode
            value = value + bmn(imode)*cos(real(ixm(imode), dp)*theta_in &
                                           - real(ixn(imode), dp)*zeta_in)
        end do
    end function expected_bmod

end program test_booz_xform_signs
