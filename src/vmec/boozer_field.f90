module boozer_field

    use constants, only: dp
    use field_base, only: field_t
    use boozer_sub, only: get_boozer_coordinates, splint_boozer_coord
    use boozer_chartmap, only: load_boozer_from_chartmap
    use boozmn_reader, only: load_boozer_from_boozmn

    implicit none
    private

    real(dp), parameter :: cm2m = 1e-2_dp
    real(dp), parameter :: gauss2tesla = 1e-4_dp

    logical, save :: initialized = .false.

    public :: boozer_field_t

    !> Field in Boozer coordinates, loadable from several file formats.
    !!
    !! Call one init_from_* loader, then fix_to_surface before any
    !! evaluation. Only one Boozer field can be active at a time (singleton —
    !! the splines live in module-global state). Calling a loader a second
    !! time aborts.
    type, extends(field_t) :: boozer_field_t
        logical :: initialized = .false.
        logical :: fixed_to_surface = .false.
        real(dp) :: fixed_stor
        real(dp) :: nfp
        real(dp) :: psi_tor_edge
        real(dp) :: R
    contains
        procedure :: init_from_vmec
        procedure :: init_from_boozmn
        procedure :: init_from_chartmap
        procedure :: evaluate
        procedure :: get_iota
        procedure :: get_covariant_components
        procedure :: fix_to_surface
        procedure :: compute_B_sqrtg_dB_dx
        procedure :: compute_B_and_dB_dx
        procedure :: compute_B_mod
        procedure :: compute_nabla_s
        procedure :: rel_accuracy_B
    end type boozer_field_t

contains

    !> Load Boozer-coordinate splines from a VMEC wout netCDF file.
    !!
    !! Reasonable defaults: radial_spline_order=5,
    !! angular_spline_order=5, grid_refinement=6.
    subroutine init_from_vmec(self, vmec_file, &
                              radial_spline_order, &
                              angular_spline_order, &
                              grid_refinement)
        use boozer_coordinates_mod, only: use_B_r
        class(boozer_field_t), intent(inout) :: self
        character(len=*), intent(in) :: vmec_file
            !! path to VMEC .nc file
        integer, intent(in), optional :: radial_spline_order
        integer, intent(in), optional :: angular_spline_order
        integer, intent(in), optional :: grid_refinement

        call assert_no_active_field("init_from_vmec")
        use_B_r = .true.
        call get_boozer_coordinates(vmec_file, &
                                    radial_spline_order, &
                                    angular_spline_order, &
                                    grid_refinement)
        call finish_init(self)

    end subroutine init_from_vmec

    !> Load Boozer-coordinate splines from a booz_xform boozmn netCDF file.
    !!
    !! The boozmn file fixes its own grid and spline orders.
    subroutine init_from_boozmn(self, boozmn_file)
        class(boozer_field_t), intent(inout) :: self
        character(len=*), intent(in) :: boozmn_file

        call assert_no_active_field("init_from_boozmn")
        call load_boozer_from_boozmn(boozmn_file)
        call finish_init(self)

    end subroutine init_from_boozmn

    !> Load Boozer-coordinate splines from an extended chartmap netCDF file.
    !!
    !! The chartmap carries its own grid and spline orders.
    subroutine init_from_chartmap(self, chartmap_file)
        class(boozer_field_t), intent(inout) :: self
        character(len=*), intent(in) :: chartmap_file

        call assert_no_active_field("init_from_chartmap")
        call load_boozer_from_chartmap(chartmap_file)
        call finish_init(self)

    end subroutine init_from_chartmap

    subroutine assert_no_active_field(caller)
        character(len=*), intent(in) :: caller

        if (initialized) error stop caller// &
            ": a Boozer field is already loaded; "// &
            "only one can be active at a time (singleton)."
    end subroutine assert_no_active_field

    subroutine finish_init(self)
        use vector_potentail_mod, only: torflux
        use new_vmec_stuff_mod, only: nper, rmajor
        class(boozer_field_t), intent(inout) :: self

        self%psi_tor_edge = -torflux*cm2m**2.0_dp*gauss2tesla
        self%nfp = real(nper, dp)
        self%R = rmajor
        initialized = .true.
        self%initialized = .true.
    end subroutine finish_init

    subroutine evaluate(self, x, bmod, sqrtg, bder, &
                        hcovar, hctrvr, hcurl)
        use vector_potentail_mod, only: torflux
        class(boozer_field_t), intent(in) :: self
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: bmod, sqrtg
        real(dp), intent(out) :: bder(3), hcovar(3), hctrvr(3), hcurl(3)

        real(dp) :: r, vartheta_B, varphi_B, &
                    A_phi, A_theta, dA_phi_dr, dA_theta_dr, &
                    d2A_phi_dr2, d3A_phi_dr3, &
                    B_vartheta_B, dB_vartheta_B, d2B_vartheta_B, &
                    B_varphi_B, dB_varphi_B, d2B_varphi_B, &
                    Bmod_B, sqrt_g_ss_B, B_r
        real(dp), dimension(3) :: dBmod_B, dB_r
        real(dp), dimension(6) :: d2Bmod_B, d2B_r

        real(dp) :: aiota, Bctrvr_theta, Bctrvr_phi, sqrtgbmod

        integer, parameter :: mode_secders = 0

        if (.not. initialized) then
            error stop "boozer_field_evaluate: field not initialized. "// &
                "Call an init_from_* loader first!"
        end if

        r = x(1)
        vartheta_B = x(2)
        varphi_B = x(3)

        call splint_boozer_coord(r, vartheta_B, varphi_B, &
                                 mode_secders, &
                                 A_theta, A_phi, dA_theta_dr, &
                                 dA_phi_dr, d2A_phi_dr2, &
                                 d3A_phi_dr3, &
                                 B_vartheta_B, dB_vartheta_B, &
                                 d2B_vartheta_B, &
                                 B_varphi_B, dB_varphi_B, &
                                 d2B_varphi_B, &
                                 Bmod_B, dBmod_B, d2Bmod_B, &
                                 B_r, dB_r, d2B_r, &
                                 sqrt_g_ss_B)

        aiota = -dA_phi_dr/dA_theta_dr

        bmod = Bmod_B
        bder = dBmod_B/Bmod_B

        sqrtg = (aiota*B_vartheta_B + B_varphi_B)/Bmod_B**2*torflux

        Bctrvr_phi = dA_theta_dr/sqrtg
        Bctrvr_theta = aiota*Bctrvr_phi
        hctrvr(1) = 0.d0
        hctrvr(2) = Bctrvr_theta/bmod
        hctrvr(3) = Bctrvr_phi/bmod

        hcovar(1) = B_r/bmod
        hcovar(2) = B_vartheta_B/bmod
        hcovar(3) = B_varphi_B/bmod

        sqrtgbmod = sqrtg*bmod
        hcurl(1) = (B_vartheta_B*bder(3) &
                    - B_varphi_B*bder(2))/sqrtgbmod
        hcurl(2) = (B_varphi_B*bder(1) - B_r*bder(3) &
                    + dB_r(3) - dB_varphi_B)/sqrtgbmod
        hcurl(3) = (B_r*bder(2) - B_vartheta_B*bder(1) &
                    + dB_vartheta_B - dB_r(2))/sqrtgbmod

        bmod = bmod*gauss2tesla
        sqrtg = sqrtg*cm2m**3.0_dp
        hcovar = hcovar*cm2m
        hctrvr = hctrvr/cm2m
        hcurl = hcurl/cm2m**2.0_dp

    end subroutine evaluate

    !> Return the rotational transform iota at flux surface stor.
    subroutine get_iota(self, stor, iota)
        use new_vmec_stuff_mod, only: nper

        class(boozer_field_t), intent(in) :: self
        real(dp), intent(in) :: stor
            !! normalized toroidal flux s, range [0, 1]
        real(dp), intent(out) :: iota

        real(dp) :: A_phi, A_theta, dA_phi_dr, dA_theta_dr
        real(dp) :: d2A_phi_dr2, d3A_phi_dr3
        real(dp) :: B_vartheta_B, B_varphi_B
        real(dp) :: dB_vartheta_B, d2B_vartheta_B
        real(dp) :: dB_varphi_B, d2B_varphi_B
        real(dp) :: Bmod_B, sqrt_g_ss_B, B_r
        real(dp), dimension(3) :: dBmod_B, dB_r
        real(dp), dimension(6) :: d2Bmod_B, d2B_r

        call splint_boozer_coord(stor, 0.0_dp, 0.0_dp, 0, &
                                 A_theta, A_phi, dA_theta_dr, &
                                 dA_phi_dr, d2A_phi_dr2, &
                                 d3A_phi_dr3, &
                                 B_vartheta_B, dB_vartheta_B, &
                                 d2B_vartheta_B, &
                                 B_varphi_B, dB_varphi_B, &
                                 d2B_varphi_B, &
                                 Bmod_B, dBmod_B, d2Bmod_B, &
                                 B_r, dB_r, d2B_r, &
                                 sqrt_g_ss_B)

        iota = -dA_phi_dr/dA_theta_dr
    end subroutine get_iota

    !> Return covariant B_theta and B_phi for the surface fixed by fix_to_surface.
    !!
    !! These are flux-surface constants in Boozer coordinates, angle-independent, SI: T*m.
    !!
    !! Requires fix_to_surface to have been called first.
    subroutine get_covariant_components(self, B_theta_covariant, B_phi_covariant)
        class(boozer_field_t), intent(in) :: self

        real(dp), intent(out) :: B_theta_covariant, B_phi_covariant

        real(dp) :: A_phi, A_theta, dA_phi_dr, dA_theta_dr
        real(dp) :: d2A_phi_dr2, d3A_phi_dr3
        real(dp) :: dB_vartheta_B, d2B_vartheta_B
        real(dp) :: dB_varphi_B, d2B_varphi_B
        real(dp) :: Bmod_B, sqrt_g_ss_B, B_r
        real(dp), dimension(3) :: dBmod_B, dB_r
        real(dp), dimension(6) :: d2Bmod_B, d2B_r

        if (.not. self%fixed_to_surface) &
            error stop "get_covariant_components: call fix_stor first"

        call splint_boozer_coord(self%fixed_stor, 0.0_dp, 0.0_dp, 0, &
                                 A_theta, A_phi, dA_theta_dr, &
                                 dA_phi_dr, d2A_phi_dr2, &
                                 d3A_phi_dr3, &
                                 B_theta_covariant, dB_vartheta_B, &
                                 d2B_vartheta_B, &
                                 B_phi_covariant, dB_varphi_B, &
                                 d2B_varphi_B, &
                                 Bmod_B, dBmod_B, d2Bmod_B, &
                                 B_r, dB_r, d2B_r, &
                                 sqrt_g_ss_B)

        B_phi_covariant = B_phi_covariant*cm2m*gauss2tesla
        B_theta_covariant = B_theta_covariant*cm2m*gauss2tesla
    end subroutine get_covariant_components

    !> Fix the field to the flux surface at normalized toroidal flux stor.
    !!
    !! Must be called before any compute_* call or get_covariant_components.
    subroutine fix_to_surface(self, stor)
        class(boozer_field_t), intent(inout) :: self
        real(dp), intent(in) :: stor
            !! normalized toroidal flux s, range [0, 1]

        self%fixed_stor = stor
        self%fixed_to_surface = .true.
    end subroutine fix_to_surface

    subroutine compute_B_sqrtg_dB_dx(self, theta, phi, B_mod, sqrtg, &
                                     dB_dx)
        class(boozer_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod, sqrtg, dB_dx(3)
        real(dp) :: dlnB_dx(3)

        real(dp) :: x(3), hcovar(3), hctrvr(3), hcurl(3)

        if (.not. self%fixed_to_surface) &
            error stop "compute_B_sqrtg_dB_dx: call fix_stor first"
        x(1) = self%fixed_stor
        x(2) = theta
        x(3) = phi
        call self%evaluate(x, B_mod, sqrtg, dlnB_dx, hcovar, &
                           hctrvr, hcurl)
        dB_dx = dlnB_dx*B_mod
    end subroutine compute_B_sqrtg_dB_dx

    subroutine compute_B_and_dB_dx(self, theta, phi, B_mod, dB_dx)
        class(boozer_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod, dB_dx(3)
        real(dp) :: dlnB_dx(3)

        real(dp) :: x(3), dummy_sqrtg, hcovar(3), hctrvr(3), hcurl(3)

        if (.not. self%fixed_to_surface) &
            error stop "compute_B_and_dB_dx: call fix_stor first"
        x(1) = self%fixed_stor
        x(2) = theta
        x(3) = phi
        call self%evaluate(x, B_mod, dummy_sqrtg, dlnB_dx, hcovar, &
                           hctrvr, hcurl)
        dB_dx = dlnB_dx*B_mod
    end subroutine compute_B_and_dB_dx

    subroutine compute_B_mod(self, theta, phi, B_mod)
        class(boozer_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: B_mod

        real(dp) :: x(3), dummy_sqrtg, dummy_dB_dx(3)
        real(dp) :: hcovar(3), hctrvr(3), hcurl(3)

        if (.not. self%fixed_to_surface) &
            error stop "compute_B_mod: call fix_stor first"
        x(1) = self%fixed_stor
        x(2) = theta
        x(3) = phi
        call self%evaluate(x, B_mod, dummy_sqrtg, dummy_dB_dx, &
                           hcovar, hctrvr, hcurl)
    end subroutine compute_B_mod

    subroutine compute_nabla_s(self, theta, phi, nabla_s)
        class(boozer_field_t), intent(in) :: self
        real(dp), intent(in) :: theta, phi
        real(dp), intent(out) :: nabla_s

        real(dp) :: dummy(32), sqrt_g_ss

        if (.not. self%fixed_to_surface) &
            error stop "compute_nabla_s: call fix_to_surface first"

        call splint_boozer_coord(self%fixed_stor, theta, phi, 0, &
                                 dummy(1), dummy(2), dummy(3), &
                                 dummy(4), dummy(5), dummy(6), &
                                 dummy(7), dummy(8), dummy(9), &
                                 dummy(10), dummy(11), dummy(12), &
                                 dummy(13), dummy(14:16), dummy(17:22), &
                                 dummy(23), dummy(24:26), dummy(27:32), &
                                 sqrt_g_ss)
        nabla_s = sqrt_g_ss/cm2m

    end subroutine compute_nabla_s

    real(dp) function rel_accuracy_B(self)
        class(boozer_field_t), intent(in) :: self

        rel_accuracy_B = 1e-9_dp
    end function rel_accuracy_B

end module boozer_field
