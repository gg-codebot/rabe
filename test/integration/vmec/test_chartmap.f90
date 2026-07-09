!> Chartmap loader test for boozer_field_t.
!>
!> Loads input/circ_chartmap.nc via init_from_chartmap and checks Bmod at
!> three grid nodes against the values stored in that file (converted from
!> Gauss to Tesla), pinning the wiring of the libneo chartmap backend into
!> the field interface.
program test_chartmap
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_field, only: boozer_field_t
    use utils, only: not_same

    implicit none

    character(len=*), parameter :: chartmap_file = 'input/circ_chartmap.nc'

    ! Grid-node coordinates and Bmod of circ_chartmap.nc at (rho, theta, zeta)
    ! indices (10,13,13), (5,7,7), (15,19,7); stor = rho**2, Bmod in Tesla.
    integer, parameter :: n_points = 3
    real(dp), parameter :: stor(n_points) = [0.22487562326869806_dp, &
                                             0.044654362880886422_dp, &
                                             0.54332416897506919_dp]
    real(dp), parameter :: theta(n_points) = [3.1415926535897931_dp, &
                                              1.5707963267948966_dp, &
                                              4.7123889803846897_dp]
    real(dp), parameter :: zeta(n_points) = [3.1415926535897931_dp, &
                                             1.5707963267948966_dp, &
                                             1.5707963267948966_dp]
    real(dp), parameter :: Bmod_ref(n_points) = [2.2254159391411612_dp, &
                                                 1.9672547079464935_dp, &
                                                 2.0602177739241334_dp]

    ! Interpolating splines reproduce grid-node values exactly (measured 0.0
    ! reldiff locally); 1e-9 leaves headroom for compiler/FMA variation.
    real(dp), parameter :: reltol = 1.0e-9_dp
    real(dp), parameter :: abstol = 1.0e-12_dp

    type(boozer_field_t) :: field
    real(dp) :: Bmod
    integer :: point
    logical :: test_failed

    call field%init_from_chartmap(chartmap_file)
    if (.not. field%initialized) error stop 'chartmap field not initialized'

    test_failed = .false.
    do point = 1, n_points
        call field%fix_to_surface(stor(point))
        call field%compute_B_mod(theta(point), zeta(point), Bmod)
        if (not_same(Bmod, Bmod_ref(point), &
                     reltol_in=reltol, abstol_in=abstol)) then
            print *, 'FAIL: chartmap Bmod mismatch at point', point
            print *, '  computed: ', Bmod
            print *, '  reference:', Bmod_ref(point)
            print *, '  reldiff:  ', abs(Bmod - Bmod_ref(point))/Bmod_ref(point)
            test_failed = .true.
        end if
    end do

    if (test_failed) error stop

end program test_chartmap
