module params
    use constants, only: dp

    implicit none

    double precision, parameter :: c = 2.9979d10
    double precision, parameter :: e_charge = 4.8032d-10
    double precision, parameter :: p_mass = 1.6726d-24
    double precision, parameter :: ev = 1.6022d-12
    real(dp) :: rlarm
    real(dp) :: dtaumin

    integer :: ntau ! number of dtaumin in dtau

    real(dp) :: relerr = 1d-10

    ! variables to evaluate at tip: z(1..5), par_inv
    integer :: nplagr, nder, npl_half

    ! Further configuration parameters
    real(dp) :: facE_al = 1d0
    integer :: ntimstep = 10000, n_e = 2
    real(dp) :: n_d = 4

    real(dp) :: v0

    ! neglect relativistic effects by large inverse relativistic temperature
    real(dp), parameter :: rmu = 1e8
    real(dp) :: ro0

contains

    subroutine params_init(dtaumin_in, rlarm_in)
        real(dp), intent(in) :: dtaumin_in
        real(dp), intent(in), optional :: rlarm_in
        real(dp) :: E_alpha
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
        dtaumin = dtaumin_in
        print *, "dtaumin = ", dtaumin
        ntau = 100
        nplagr = 6
        nder = 0
        npl_half = nplagr/2

    end subroutine params_init

end module params
