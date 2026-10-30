program electric_rabe
    use constants, only: dp, pi
    use utils, only: linspace
    use boozer_field, only: boozer_field_t
    use fieldline_mod, only: fieldline_t, flock_of_fieldlines_t
    use make_fieldline, only: make_flock_of_fieldlines

    use precession, only: get_smallest_maximum, get_biggest_minimum
    use precession, only: set_fieldline_minima
    use grid_mod, only: integration_grid_t
    use grid_mod, only: set_integration_grids, fieldline_for_precession_t
    use grid_mod, only: compute_bounce_integrals
    use field_instance, only: initialize_field_instance
    use splines_instance, only: initialize_splines
    use splines_instance, only: initialize_prefactor
    use splines_instance, only: initialize_startup
    use splines_instance, only: initialize_radial_drift_spline
    use splines_instance, only: get_flux_mode
    use fourier, only: real_ft
    use fieldline_labels, only: modes_t

    use surface_average_mod, only: surface_average_t, calc_surface_averages
    use coefficients, only: calc_nu_star_crit
    use coefficients, only: calc_finite_boundary_layer_correction
    use netcdf_mod, only: netcdf_t
    use fieldline_labels, only: fieldline_modes_t, fourier_transform_over_label
    use git_version, only: git_hash

    use read_file, only: read_namelist
    use read_file, only: field_file, &
                         M_pol, &
                         N_tor, &
                         s_tor, &
                         sign_sqrtg, &
                         max_n_fieldlines, &
                         should_calc_shaing_callen, &
                         n_eta
    use read_precession, only: read_precession_namelist
    use read_precession, only: Omega_hat, &
                               nu_star

    implicit none

    character(len=*), parameter :: input_file = "rabe.in"
    character(len=*), parameter :: output_file = "electric_rabe.nc"

    type(boozer_field_t) :: field

    integer :: n_stor
    integer :: this

    real(dp) :: R ![m]
    type(surface_average_t) :: average
    real(dp) :: dr_dAtheta ![rad/Tm/]
    real(dp) :: iota
    real(dp) :: B_theta_covariant, B_phi_covariant
    real(dp) :: nfp

    integer :: n_fieldlines
    type(flock_of_fieldlines_t) :: flock
    type(fieldline_t), dimension(:), allocatable :: fieldlines
    logical :: too_strong_violation

    real(dp) :: covariant_factor
    real(dp), dimension(:), allocatable :: nu_star_crit
    integer, dimension(:), allocatable :: err_flag
    real(dp) :: helical_factor

    type(fieldline_modes_t) :: modes
    type(fieldline_for_precession_t), dimension(:), allocatable :: fieldlines_precession
    real(dp) :: phi_bottom, B_bottom, lowest_B_max, eta_t, eta_c
    type(integration_grid_t) :: grid

    real(dp), dimension(:, :), allocatable :: radial_drift_weighted
    real(dp), dimension(:, :), allocatable :: radial_drift_cos
    real(dp), dimension(:, :), allocatable :: radial_drift_sin
    real(dp), dimension(:), allocatable :: bounce_time_weighted
    real(dp), dimension(:), allocatable :: I_j
    real(dp), dimension(:), allocatable :: poloidal_drift_weighted
    real(dp), dimension(:), allocatable :: electric_drift_weighted
    real(dp), dimension(:), allocatable :: magnetic_drift_weighted
    logical, parameter :: ignore_magnetic_drift = .true.

    real(dp), dimension(:), allocatable :: flux_mode
    real(dp) :: mode_factor
    type(modes_t) :: g_off_modes

    integer :: n_modes
    integer :: idx, i_grid
    integer :: file_unit
    character(len=*), parameter :: save_dir = "saved_bounce_integrals/"
    character(len=256) :: filename
    logical :: file_exists
    character(len=256) :: header_line
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
        fieldlines = flock%fieldlines

        covariant_factor = (B_phi_covariant + B_theta_covariant*iota)
        call calc_surface_averages(flock, average)
        dr_dAtheta = sign_sqrtg/(average%nabla_s*field%psi_tor_edge)
        R = field%R

! --------------------------------------------------------------------------
        allocate (fieldlines_precession(n_fieldlines))
        fieldlines_precession%fieldline_t = fieldlines
        fieldlines_precession%M_pol = M_pol
        fieldlines_precession%N_tor = N_tor
        fieldlines_precession%nfp = nfp
        eta_c = 1.0_dp/get_smallest_maximum(fieldlines_precession)
        call set_fieldline_minima(field, fieldlines_precession)
        eta_t = 1.0_dp/get_biggest_minimum(fieldlines_precession)
        call set_integration_grids(eta_t, eta_c, grid)
        fieldlines_precession%grid = grid

        call initialize_field_instance(field)
        allocate (radial_drift_weighted(n_fieldlines, grid%n))
        allocate (bounce_time_weighted(grid%n))
        allocate (I_j(grid%n))
        allocate (magnetic_drift_weighted(grid%n))
        allocate (electric_drift_weighted(grid%n))
        bounce_time_weighted = 0.0_dp
        I_j = 0.0_dp
        magnetic_drift_weighted = 0.0_dp

        do idx = 1, n_fieldlines
            write (filename, '(A,A,I0,A)') save_dir, 'bounce_fieldline_', idx, '.dat'
            inquire (file=trim(filename), exist=file_exists)

            if (file_exists) then
                print *, "Reading bounce integrals from ", trim(filename)
                open (newunit=file_unit, file=trim(filename), &
                      status='old', action='read')
                read (file_unit, '(A)') header_line
                do i_grid = 1, fieldlines_precession(idx)%grid%n
                    read (file_unit, *) &
                        fieldlines_precession(idx)%grid%t(i_grid), &
                        fieldlines_precession(idx)%grid%eta(i_grid), &
                        fieldlines_precession(idx)%grid%bounce_time_weighted(i_grid), &
                        fieldlines_precession(idx)%grid%I_j(i_grid), &
                        fieldlines_precession(idx)%grid%radial_drift_weighted(i_grid), &
                        fieldlines_precession(idx)%grid%poloidal_drift_weighted(i_grid)
                end do
                close (file_unit)
            else
                call compute_bounce_integrals(field, &
                                              fieldlines_precession(idx), &
                                              s_tor(this), &
                                              fieldlines_precession(idx)%grid)
                open (newunit=file_unit, file=trim(filename), &
                      status='replace')
                write (file_unit, '(A)') &
                    '# t  eta  bounce_time_weighted' &
                    //'  I_j  radial_drift_weighted' &
                    //'  poloidal_drift_weighted'
                do i_grid = 1, fieldlines_precession(idx)%grid%n
                    write (file_unit, '(6ES20.12)') &
                        fieldlines_precession(idx)%grid%t(i_grid), &
                        fieldlines_precession(idx)%grid%eta(i_grid), &
                        fieldlines_precession(idx)%grid%bounce_time_weighted(i_grid), &
                        fieldlines_precession(idx)%grid%I_j(i_grid), &
                        fieldlines_precession(idx)%grid%radial_drift_weighted(i_grid), &
                        fieldlines_precession(idx)%grid%poloidal_drift_weighted(i_grid)
                end do
                close (file_unit)
            end if

   radial_drift_weighted(idx, :) = fieldlines_precession(idx)%grid%radial_drift_weighted
            bounce_time_weighted = bounce_time_weighted + &
                                   fieldlines_precession(idx)%grid%bounce_time_weighted
            I_j = I_j + fieldlines_precession(idx)%grid%I_j
            magnetic_drift_weighted = magnetic_drift_weighted + &
                                 fieldlines_precession(idx)%grid%poloidal_drift_weighted
        end do
        if (ignore_magnetic_drift) then
            magnetic_drift_weighted = 0.0_dp
        else
            magnetic_drift_weighted = magnetic_drift_weighted/n_fieldlines
        end if
        bounce_time_weighted = bounce_time_weighted/n_fieldlines
        I_j = I_j/n_fieldlines
        allocate (poloidal_drift_weighted(grid%n))

        n_modes = n_fieldlines/2 + 1
        allocate (radial_drift_cos(n_modes, grid%n))
        allocate (radial_drift_sin(n_modes, grid%n))
        do idx = 1, grid%n
            call real_ft(fieldlines_precession%xi_0, &
                         radial_drift_weighted(:, idx), &
                         radial_drift_cos(:, idx), &
                         radial_drift_sin(:, idx))
        end do

        allocate (flux_mode(n_modes))
        allocate (lambda_off_modes(n_modes, n_nu, n_Omega_hat))

        do idx_Omega = 1, n_Omega_hat
            print *, "idx_Omega: ", idx_Omega, " / ", n_Omega_hat
            electric_drift_weighted = -Omega_hat(idx_Omega)*bounce_time_weighted
            poloidal_drift_weighted = magnetic_drift_weighted + electric_drift_weighted

            call initialize_splines(grid%t, &
                                    grid%eta, &
                                    I_j, &
                                    poloidal_drift_weighted)

            do id_nu = 1, n_nu
                l_c = pi*R/(2.0_dp*nu_star(id_nu))
                flux_mode(1) = 0.0_dp
                do idx = 2, n_modes
                    mode_factor = 0.5_dp*real(idx - 1, dp)*l_c*nfp/(M_pol*iota - N_tor)
                    call initialize_prefactor(mode_factor)
                    call initialize_startup(grid%t, grid%eta, I_j, mode_factor)
                   call initialize_radial_drift_spline(grid%t, radial_drift_sin(idx, :))
                    call get_flux_mode(grid%t(1), grid%t(grid%n), flux_mode(idx))
                end do
                call get_g_modes_from_fieldlines(flock, l_c, g_off_modes, &
                                                 covariant_factor)
        lambda_off_modes(:, id_nu, idx_Omega) = pi*g_off_modes%sin_coeffs*flux_mode/average%normalization
               lambda_off(id_nu, idx_Omega) = sum(lambda_off_modes(:, id_nu, idx_Omega))

            end do
        end do
        lambda_off = lambda_off*0.75_dp/covariant_factor*sign_sqrtg/average%nabla_s
 lambda_off_modes = lambda_off_modes*0.75_dp/covariant_factor*sign_sqrtg/average%nabla_s

        deallocate (flux_mode)
        deallocate (fieldlines_precession)
        deallocate (radial_drift_weighted)
        deallocate (bounce_time_weighted)
        deallocate (I_j)
        deallocate (poloidal_drift_weighted)
        deallocate (electric_drift_weighted)
        deallocate (magnetic_drift_weighted)
        deallocate (radial_drift_cos)
        deallocate (radial_drift_sin)

! --------------------------------------------------------------------------

        nu_star_crit(this) = calc_nu_star_crit(flock, R)

        print *, "s_tor: ", s_tor(this)

        if (allocated(fieldlines)) deallocate (fieldlines)

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
        real(dp), dimension(:), allocatable :: theta_mid, g_off

        character(len=1024) :: label

        call fourier_transform_over_label(flock, modes)
        call calc_surface_averages(flock, averages)

        max_mode = size(modes%delta_eta%cos_coeffs, dim=1)
        if (.not. allocated(g_off_modes%sin_coeffs)) then
            call allocate_modes(g_off_modes, max_mode)
        end if

        prefactor_A = 2.0_dp*sqrt(covariant_factor*flock%eta_b* &
                                  flock%I_ref/l_c)
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
