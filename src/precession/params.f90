module params
    use, intrinsic :: iso_fortran_env, only: int8
    use constants, only: dp, pi

    implicit none

    double precision, parameter :: c = 2.9979d10
    double precision, parameter :: e_charge = 4.8032d-10
    double precision, parameter :: e_mass = 9.1094d-28
    double precision, parameter :: p_mass = 1.6726d-24
    double precision, parameter :: ev = 1.6022d-12
    integer :: nper = 1000, ntestpart = 1024
    integer :: zstart_dim1 = 5
    real(dp) :: dphi, phibeg = 0d0, bmod00, rlarm, bmax, bmin
    real(dp) :: tau, dtau, dtaumin, xi
    real(dp) :: RT0, R0i, cbfi, bz0i, bf0, rbig
    real(dp), dimension(1024) :: sbeg = 0.5d0
    real(dp) :: thetabeg = 0.0d0
    real(dp), dimension(:), allocatable :: bstart, volstart
    !where npoi is used, add a dimension at the end for sbeg
    real(dp), dimension(:, :), allocatable :: xstart
    real(dp), dimension(:, :), allocatable :: zstart, zend
    real(dp), dimension(:), allocatable :: confpart_trap, confpart_pass
    real(dp), dimension(:), allocatable :: times_lost
    real(dp) :: contr_pp = -1d0
    integer :: ibins
    logical :: generate_start_only = .false.
    integer :: startmode = 1
    real(dp) :: grid_density = 0d0
    logical :: special_ants_file = .false.

    integer :: ntau ! number of dtaumin in dtau
    integer(8) :: n_microsteps_total = 0_8
    character(16) :: macrostep_time_grid = 'linear'

    integer :: kpart = 0 ! progress counter for particles

    real(dp) :: relerr = 1d-10

    real(dp), allocatable :: trap_par(:), perp_inv(:)
    integer, allocatable :: iclass(:, :)
    logical, allocatable :: class_passing(:), class_lost(:)

    ! variables to evaluate at tip: z(1..5), par_inv
    integer :: nplagr, nder, npl_half
    integer :: norbper, nfp
    real(dp) :: fper

    real(dp) :: tcut = -1d0
    integer :: ntcut
    integer :: nturns = 8
    logical :: class_plot = .false.    !<=AAA
    real(dp) :: cut_in_per = 0d0        !<=AAA

    logical :: fast_class = .false.
    !if .true. quit immeadiately after fast classification

    ! Colliding with D-T reactor plasma. TODO: Make configurable
    logical :: swcoll = .false.
    real(dp) :: am1 = 2.0d0, am2 = 3.0d0, Z1 = 1.0d0, Z2 = 1.0d0, &
                densi1 = 0.5d14, densi2 = 0.5d14, tempi1 = 1.0d4, tempi2 = 1.0d4, &
                tempe = 1.0d4
    real(dp) :: dchichi, slowrate, dchichi_norm, slowrate_norm
    logical :: deterministic = .false.

    ! Further configuration parameters
    integer :: notrace_passing = 0
    real(dp) :: facE_al = 1d0, trace_time = 2d-3
    integer :: ntimstep = 10000, npoiper = 100, npoiper2 = 128, n_e = 2
    real(dp) :: n_d = 4

    real(dp) :: v0

    logical :: debug = .false.
    logical :: output_results_netcdf = .false.
    integer :: ierr

    integer :: batch_size = 2000000000  ! Initialize large so batch mode is not default
    integer :: ran_seed = 12345
    integer :: num_surf = 1
    logical :: reuse_batch = .false.
    integer, dimension(:), allocatable :: idx
    real(dp), parameter :: rmu = 1e8
    real(dp) :: ro0

    character(1000) :: field_input = ''
    character(1000) :: coord_input = ''
    character(1000) :: wall_input = ''
    character(16) :: wall_units = 'm'
    integer :: integ_coords = -1000  ! Sentinel: -1000 means user did not set it

contains

    subroutine params_init(nfperiods, dtaumin_in, rlarm_in)
        real(dp), intent(in) :: nfperiods, dtaumin_in
        real(dp), intent(in), optional :: rlarm_in
        real(dp) :: E_alpha
        integer :: L1i
        real(dp) :: weight_sum, cumul_weight, w
        integer :: i, nintv
        integer(8) :: kt_target, kt_prev
        real(dp) :: m2cm = 1e2, tesla2gauss = 1e4

        E_alpha = 3.5d6/facE_al
        ! set alpha energy, velocity, and Larmor radius
        v0 = sqrt(2.d0*E_alpha*ev/(n_d*p_mass))
        if (present(rlarm_in)) then
            rlarm = rlarm_in*m2cm*tesla2gauss
        else
            rlarm = v0*n_d*p_mass*c/(n_e*e_charge)
        end if
        ro0 = rlarm
        ! Neglect relativistic effects by large inverse relativistic temperature
        ! normalized slowing down time:
        tau = trace_time*v0
        ! normalized time step:
        dtau = tau/dble(ntimstep - 1)
        ! parameters for the vacuum chamber:
        L1i = int(nfperiods)
        ! rt0 = rmajor*1e2
        ! rbig = rt0
        ! field line integration step step over phi (to check chamber wall crossing)
        dphi = 2.d0*pi/(L1i*npoiper)
        ! orbit integration time step (to check chamber wall crossing)
        ! dtaumin = 2.d0*pi*rbig/npoiper2
        dtaumin = dtaumin_in
        print *, "dtaumin = ", dtaumin
        ntau = 100
        fper = 2d0*pi/dble(L1i)
        nplagr = 6
        nder = 0
        npl_half = nplagr/2

    end subroutine params_init

end module params
