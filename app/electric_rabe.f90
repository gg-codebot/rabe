program electric_rabe
    use constants, only: dp, pi
    use boozer_field, only: boozer_field_t
    use fieldline_mod, only: flock_of_fieldlines_t
    use make_fieldline, only: make_flock_of_fieldlines

    use precession, only: bounce_integrals_t
    use precession, only: make_fieldlines_for_precession
    use precession, only: assemble_bounce_integrals
    use precession, only: compute_flux_modes
    use grid_mod, only: fieldline_for_precession_t
    use grid_mod, only: compute_bounce_integrals
    use saved_bounce_integrals, only: bounce_integrals_header_t
    use saved_bounce_integrals, only: write_bounce_integrals_header
    use saved_bounce_integrals, only: append_bounce_integrals, read_bounce_integrals
    use field_instance, only: initialize_field_instance
    use fieldline_labels, only: modes_t

    use surface_average_mod, only: surface_average_t, calc_surface_averages
    use coefficients, only: calc_nu_star_crit
    use netcdf_mod, only: netcdf_t
    use fieldline_labels, only: fieldline_modes_t, fourier_transform_over_label
    use git_version, only: git_hash

    use read_file, only: read_namelist
    use read_file, only: field_file, &
                         M_pol, &
                         N_tor, &
                         s_tor, &
                         sign_sqrtg, &
                         max_n_fieldlines
    use read_precession, only: read_precession_namelist
    use read_precession, only: Omega_hat, &
                               nu_star, &
                               n_eta

    implicit none

    character(len=*), parameter :: input_file = "rabe.in"
    character(len=*), parameter :: output_file = "electric_rabe.nc"
    character(len=*), parameter :: bounce_integrals_file = "saved_bounce_integrals.dat"

    type(boozer_field_t) :: field

    integer :: n_stor
    integer :: this

    real(dp) :: R ![m]
    type(surface_average_t) :: average
    real(dp) :: iota
    real(dp) :: B_theta_covariant, B_phi_covariant
    real(dp) :: nfp

    integer :: n_fieldlines
    type(flock_of_fieldlines_t) :: flock

    real(dp) :: covariant_factor
    real(dp), dimension(:), allocatable :: nu_star_crit
    integer, dimension(:), allocatable :: err_flag

    type(fieldline_for_precession_t), dimension(:), allocatable :: fieldlines_precession
    type(bounce_integrals_t) :: bounce_integrals

    real(dp), dimension(:), allocatable :: flux_mode
    type(modes_t) :: g_off_modes

    integer :: n_modes
    integer :: idx
    logical :: file_exists
    type(bounce_integrals_header_t) :: bounce_integrals_header
    logical, dimension(:), allocatable :: is_complete
    integer :: n_Omega_hat
    integer :: idx_Omega
    integer :: n_nu
    integer :: id_nu
    real(dp) :: l_c
    real(dp), dimension(:, :), allocatable :: lambda_off
    real(dp), dimension(:, :, :), allocatable :: lambda_off_modes

    type(netcdf_t) :: nc_output
    character(len=*), parameter :: dim_name = "surface"
    character(len=1024) :: description

    call read_namelist(input_file)
    call read_precession_namelist(input_file)

    n_stor = size(s_tor)
    if (n_stor > 1) then
        error stop "Error: only one surface allowed!"
    end if
    allocate (nu_star_crit(n_stor))
    allocate (err_flag(n_stor))
    n_Omega_hat = size(Omega_hat)
    n_nu = size(nu_star)
    allocate (lambda_off(n_nu, n_Omega_hat))

    call field%init_from_vmec(field_file, grid_refinement=6)
    do this = 1, n_stor
        call field%fix_to_surface(s_tor(this))
        call field%get_iota(s_tor(this), iota)
        call field%get_covariant_components(B_theta_covariant, &
                                            B_phi_covariant)
        nfp = field%nfp

        call make_flock_of_fieldlines(flock, max_n_fieldlines, iota, &
                                      field, M_pol, N_tor, nfp, &
                                      err_flag(this))
        iota = flock%iota
        n_fieldlines = size(flock%fieldlines)

        covariant_factor = (B_phi_covariant + B_theta_covariant*iota)
        call calc_surface_averages(flock, average)
        R = field%R

! --------------------------------------------------------------------------
        call make_fieldlines_for_precession(field, flock, fieldlines_precession, n_eta)
        call initialize_field_instance(field)

        bounce_integrals_header = &
            bounce_integrals_header_t(git_hash=git_hash, &
                                      field_file=trim(field_file), &
                                      s_tor=s_tor(this), &
                                      M_pol=M_pol, &
                                      N_tor=N_tor, &
                                      n_grid=n_eta, &
                                      max_n_fieldlines=max_n_fieldlines, &
                                      n_fieldlines=n_fieldlines, &
                                      eta_c=fieldlines_precession(1)%grid%eta(n_eta), &
                                      eta_t=fieldlines_precession(1)%grid%eta(1))
        allocate (is_complete(n_fieldlines))
        inquire (file=bounce_integrals_file, exist=file_exists)
        if (file_exists) then
            call read_bounce_integrals(bounce_integrals_file, &
                                       bounce_integrals_header, &
                                       fieldlines_precession, is_complete)
            print *, "Read ", count(is_complete), " of ", n_fieldlines, &
                " fieldlines from ", bounce_integrals_file
        else
            call write_bounce_integrals_header(bounce_integrals_file, &
                                               bounce_integrals_header)
            is_complete = .false.
        end if

        do idx = 1, n_fieldlines
            if (.not. is_complete(idx)) then
                call compute_bounce_integrals(field, &
                                              fieldlines_precession(idx), &
                                              s_tor(this), &
                                              fieldlines_precession(idx)%grid)
                call append_bounce_integrals(bounce_integrals_file, idx, &
                                             fieldlines_precession(idx)%grid)
            end if
        end do
        call assemble_bounce_integrals(fieldlines_precession, bounce_integrals)
        n_modes = bounce_integrals%n_modes

        allocate (lambda_off_modes(n_modes, n_nu, n_Omega_hat))

        do idx_Omega = 1, n_Omega_hat
            print *, "idx_Omega: ", idx_Omega, " / ", n_Omega_hat
            do id_nu = 1, n_nu
                l_c = pi*R/(2.0_dp*nu_star(id_nu))
                call compute_flux_modes(bounce_integrals, l_c, Omega_hat(idx_Omega), &
                                        flux_mode)
                call get_g_modes_from_fieldlines(flock, l_c, g_off_modes, &
                                                 covariant_factor)
        lambda_off_modes(:, id_nu, idx_Omega) = pi*g_off_modes%sin_coeffs*flux_mode/average%normalization
               lambda_off(id_nu, idx_Omega) = sum(lambda_off_modes(:, id_nu, idx_Omega))

            end do
        end do
        lambda_off = lambda_off*0.75_dp/covariant_factor*sign_sqrtg/average%nabla_s
 lambda_off_modes = lambda_off_modes*0.75_dp/covariant_factor*sign_sqrtg/average%nabla_s

        deallocate (fieldlines_precession)
        deallocate (is_complete)

! --------------------------------------------------------------------------

        nu_star_crit(this) = calc_nu_star_crit(flock, R)

        print *, "s_tor: ", s_tor(this)

    end do

    call nc_output%create(output_file)
    call nc_output%add_global_attribute("title", &
                                        "asymptotic bootstrap coefficient lambda_bB")
    call nc_output%add_global_attribute("definition", &
                                        "lambda^{off}_bB = ?")
    call nc_output%add_global_attribute("git_hash", git_hash)
    call nc_output%def_dim(dim_name, n_stor)
    call nc_output%add_real_1d("s_tor", dim_name)
    call nc_output%add_attr("s_tor", "long_name", &
                            "normalized toroidal flux label")
    call nc_output%add_attr("s_tor", "unit", &
                            "[1]")

    call nc_output%def_dim("collisionality", n_nu)
    call nc_output%add_real_1d("nu_star", "collisionality")
    call nc_output%def_dim("precession_frequency", n_Omega_hat)
    call nc_output%add_real_1d("Omega_hat", "precession_frequency")

    call nc_output%add_real_2d("lambda_off", "collisionality", "precession_frequency")
    call nc_output%add_attr("lambda_off", "long_name", &
                            "geometrical factor")
    call nc_output%add_attr("lambda_off", "unit", &
                            "[1]")

    call nc_output%def_dim("label_mode", n_modes)
    call nc_output%add_real_1d("mode_numbers", "label_mode")
    call nc_output%add_attr("mode_numbers", "long_name", &
                            "Fourier mode numbers")
    call nc_output%add_real_3d("lambda_off_modes", &
                               "label_mode", "collisionality", &
                               "precession_frequency")
    call nc_output%add_attr("lambda_off_modes", "long_name", &
                            "per-mode contribution to lambda_off")
    call nc_output%add_attr("lambda_off_modes", "unit", &
                            "[1]")

    call nc_output%add_real_1d("nu_star_crit", dim_name)
    call nc_output%add_attr("nu_star_crit", "long_name", &
                          "lower collisionality limit for validity of asymptotic model")
    call nc_output%add_attr("nu_star_crit", "unit", &
                            "[1]")
    call nc_output%add_int_1d("err_flag", dim_name)
    call nc_output%add_attr("err_flag", "long_name", &
                            "1 if violation of omnigeneity is too strong, 0 otherwise")
    call nc_output%add_real("R")
    call nc_output%add_attr("R", "long_name", &
                            "major radius")
    call nc_output%add_attr("R", "unit", "[m]")
    write (description, "(A,A,A)") "defines the reference length scale to convert ", &
        "to dimensionless quantities i.e. defines coeffients in respect to ", &
        "nu_star = pi*R/(2*mean_free_path) = pi*R*deflection_frequency/particle_speed"
    call nc_output%add_attr("R", "definition", description)
    call nc_output%write_real_2d("lambda_off", lambda_off)
    call nc_output%write_real_1d("mode_numbers", g_off_modes%mode_numbers)
    call nc_output%write_real_3d("lambda_off_modes", lambda_off_modes)
    call nc_output%write_real_1d("Omega_hat", Omega_hat)
    call nc_output%write_real_1d("nu_star", nu_star)
    call nc_output%write_real_1d("nu_star_crit", nu_star_crit)
    call nc_output%write_int_1d("err_flag", err_flag)
    call nc_output%write_real_1d("s_tor", s_tor)
    call nc_output%write_real("R", field%R)
    call nc_output%close()

    if (allocated(lambda_off_modes)) deallocate (lambda_off_modes)
    if (allocated(nu_star_crit)) deallocate (nu_star_crit)
    if (allocated(err_flag)) deallocate (err_flag)

contains

    subroutine get_g_modes_from_fieldlines(flock, l_c, g_off_modes, covariant_factor)
        use fit_functions, only: S_A, S_B
        use fieldline_labels, only: allocate_modes
        type(flock_of_fieldlines_t), intent(in) :: flock
        real(dp), intent(in) :: l_c
        real(dp), intent(in) :: covariant_factor
        type(modes_t), intent(out) :: g_off_modes

        type(fieldline_modes_t) :: modes
        type(surface_average_t) :: averages

        real(dp) :: prefactor_A

        integer :: max_mode

        call fourier_transform_over_label(flock, modes)
        call calc_surface_averages(flock, averages)

        max_mode = size(modes%delta_eta%cos_coeffs, dim=1)
        if (.not. allocated(g_off_modes%sin_coeffs)) then
            call allocate_modes(g_off_modes, max_mode)
        end if

        ! flock%I_ref already contains the covariant factor (857dc93)
        prefactor_A = 2.0_dp*sqrt(flock%eta_b*flock%I_ref/l_c)
        g_off_modes%sin_coeffs = prefactor_A*modes%delta_aspect_ratio%cos_coeffs* &
                                 S_A(flock%iota_p*modes%delta_aspect_ratio%mode_numbers)
        g_off_modes%sin_coeffs = g_off_modes%sin_coeffs + &
                                 modes%delta_eta%cos_coeffs* &
                                 S_B(flock%iota_p*modes%delta_eta%mode_numbers)
        g_off_modes%sin_coeffs = g_off_modes%sin_coeffs* &
                                 averages%B_squared/averages%lambda_b* &
                                 l_c*0.5_dp
    end subroutine get_g_modes_from_fieldlines
end program electric_rabe
