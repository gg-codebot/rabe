!> Test that the field implementation matches the SIMPLE code for a few points in the plasma.
program test_against_simple
    use constants, only: dp
    use boozer_field, only: boozer_field_t
    use utils, only: not_same

    implicit none

    !> Results of splines used depens on the compile options, so a higher tolerance is used here.
    real(dp) :: reltol = 1e-6, abstol = 1e-11
    character(len=*), parameter :: nc_filename = "input/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"

    ! as SIMPLE uses CGS units, but we use SI
    real(dp), parameter :: cm2m = 1e-2_dp
    real(dp), parameter :: gauss2tesla = 1e-4_dp

    type(boozer_field_t) :: field

    integer, parameter :: n_cases = 5
    real(dp) :: stor(n_cases), theta(n_cases), phi(n_cases)
    real(dp) :: x(3)
    real(dp) :: bmod, sqrtg
    real(dp) :: bder(3), hcovar(3)
    real(dp) :: hctrvr(3), hcurl(3)
    real(dp) :: bmod_ref(n_cases), sqrtg_ref(n_cases)
    real(dp) :: bder_ref(n_cases, 3), hcovar_ref(n_cases, 3)
    real(dp) :: hctrvr_ref(n_cases, 3), hcurl_ref(n_cases, 3)
    integer :: case

    logical :: test_failed

    call field%init_from_vmec(nc_filename, &
                              radial_spline_order=5, &
                              angular_spline_order=5, &
                              grid_refinement=3)
    test_failed = .false.
    stor = [0.1_dp, 0.3_dp, 0.5_dp, 0.7_dp, 0.9_dp]
    theta = [0.0_dp, 1.0_dp, 3.14_dp, 0.5_dp, 2.0_dp]
    phi = [0.0_dp, 0.5_dp, 1.57_dp, 2.0_dp, 3.0_dp]

    bmod_ref = [5.7024014488290690e+04_dp, &
                5.7010408592594038e+04_dp, &
                6.3742758213874906e+04_dp, &
                5.4218536572575875e+04_dp, &
                6.1083220187449493e+04_dp]

    sqrtg_ref = [-1.5711629692175472e+07_dp, &
                 -1.5719091463432141e+07_dp, &
                 -1.2574007504120229e+07_dp, &
                 -1.7379594198471960e+07_dp, &
                 -1.3692768685651677e+07_dp]

    bder_ref(1,:) = [-1.7906807464749691E-01_dp, -7.4166439213088501E-14_dp, -1.0958457393698770E-13_dp]
    bder_ref(2,:) = [-6.5375422685029619E-02_dp,  4.8901137158983657E-02_dp, -4.2593419771991015E-04_dp]
    bder_ref(3,:) = [ 7.3646069743459427E-02_dp,  1.3567156911686458E-04_dp, -2.4743000679278641E-05_dp]
    bder_ref(4,:) = [-6.5308328098055424E-02_dp,  3.9538079263501176E-02_dp, -6.8133243699882112E-04_dp]
    bder_ref(5,:) = [ 1.2755344859923703E-02_dp,  1.0110736983359410E-01_dp, -7.9541221413346563E-04_dp]

    !> reference values field evaluation ignoring B_r
    ! hcovar_ref(1,:) = [ 0.0000000000000000E+00_dp, -9.4970349047515505E-04_dp,  1.1151159177109234E+03_dp]
    ! hcovar_ref(2,:) = [ 0.0000000000000000E+00_dp, -3.8520822192834354E-04_dp,  1.1153790773146934E+03_dp]
    ! hcovar_ref(3,:) = [ 0.0000000000000000E+00_dp, -1.8842726486691493E-04_dp,  9.9757472565834667E+02_dp]
    ! hcovar_ref(4,:) = [ 0.0000000000000000E+00_dp, -9.8956245690572466E-05_dp,  1.1728116578017384E+03_dp]
    ! hcovar_ref(5,:) = [ 0.0000000000000000E+00_dp, -1.0129726085459463E-05_dp,  1.0410079194890063E+03_dp]

    hcovar_ref(1,:) = [-9.4330617146665213E-13_dp, -9.4970349047515505E-04_dp,  1.1151159177109234E+03_dp]
    hcovar_ref(2,:) = [-4.9526222447404823E-02_dp, -3.8520822192834354E-04_dp,  1.1153790773146934E+03_dp]
    hcovar_ref(3,:) = [ 8.2784211500102169E-04_dp, -1.8842726486691493E-04_dp,  9.9757472565834667E+02_dp]
    hcovar_ref(4,:) = [-5.9028395703390649E-02_dp, -9.8956245690572466E-05_dp,  1.1728116578017384E+03_dp]
    hcovar_ref(5,:) = [ 2.8027264710120063E-02_dp, -1.0129726085459463E-05_dp,  1.0410079194890063E+03_dp]

    hctrvr_ref(1,:) = [0.0000000000000000E+00_dp,  3.7881309780855648E-04_dp,  8.9676807933375402E-04_dp]
    hctrvr_ref(2,:) = [0.0000000000000000E+00_dp,  3.7723748043895933E-04_dp,  8.9655630597133623E-04_dp]
    hctrvr_ref(3,:) = [0.0000000000000000E+00_dp,  4.2022807275813337E-04_dp,  1.0024312499722556E-03_dp]
    hctrvr_ref(4,:) = [0.0000000000000000E+00_dp,  3.5621175223770338E-04_dp,  8.5265185470933104E-04_dp]
    hctrvr_ref(5,:) = [0.0000000000000000E+00_dp,  4.0002199335749606E-04_dp,  9.6060748946365118E-04_dp]

    !> reference values field evaluation ignoring B_r
    ! hcurl_ref(1,:) = [-5.2638894003780156E-18_dp,  1.2709891638015171E-05_dp,  1.5769403006971614E-10_dp]
    ! hcurl_ref(2,:) = [ 3.4698764369963150E-06_dp,  4.6384590776612628E-06_dp, -9.1597067639039990E-11_dp]
    ! hcurl_ref(3,:) = [ 1.0763674479658840E-08_dp, -5.8430548976060074E-06_dp, -5.2040711293125020E-11_dp]
    ! hcurl_ref(4,:) = [ 2.6681129427048551E-06_dp,  4.4072755802204980E-06_dp, -3.2761515796050763E-11_dp]
    ! hcurl_ref(5,:) = [ 7.6867998812914078E-06_dp, -9.6986625274201992E-07_dp, -4.5236550526892759E-11_dp]

    hcurl_ref(1,:) = [-5.2638894003780156E-18_dp,  1.2740530997096499E-05_dp, -8.6197541867643010E-09_dp]
    hcurl_ref(2,:) = [ 3.4698764369963150E-06_dp,  4.6481422398506964E-06_dp, -1.7483693891917036E-08_dp]
    hcurl_ref(3,:) = [ 1.0763674479658840E-08_dp, -5.7798304929806501E-06_dp, -9.7743856721982719E-09_dp]
    hcurl_ref(4,:) = [ 2.6681129427048551E-06_dp,  4.3919431197038171E-06_dp, -5.0485824545366451E-09_dp]
    hcurl_ref(5,:) = [ 7.6867998812914078E-06_dp, -9.5911521268646427E-07_dp,  7.0764342860889441E-09_dp]

    do case = 1, n_cases
        x = [stor(case), theta(case), phi(case)]
        call field%evaluate(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        bmod = bmod/gauss2tesla
        if (not_same(bmod, bmod_ref(case), reltol_in=reltol, abstol_in=abstol)) then
            print *, "---------------------------------------------------"
            print *, "test_against_simple failed: B for case ", case
            print *, "B: ", bmod
            print *, "SIMPLE: ", bmod_ref(case)
            print *, "Relative error: ", abs(bmod - bmod_ref(case))/abs(bmod_ref(case))
            print *, "Absolute error: ", abs(bmod - bmod_ref(case))
            test_failed = .true.
        end if
        sqrtg = sqrtg/(cm2m**3.0_dp)
        if (not_same(sqrtg, sqrtg_ref(case), reltol_in=reltol, abstol_in=abstol)) then
            print *, "---------------------------------------------------"
            print *, "test_against_simple failed: sqrt(g) for case ", case
            print *, "sqrt(g): ", sqrtg
            print *, "SIMPLE: ", sqrtg_ref(case)
            print *, "Relative error: ", abs(sqrtg - sqrtg_ref(case)) &
                /abs(sqrtg_ref(case))
            print *, "Absolute error: ", abs(sqrtg - sqrtg_ref(case))
            test_failed = .true.
        end if
        if (not_same(bder, bder_ref(case, :), reltol_in=reltol, abstol_in=abstol)) then
            print *, "---------------------------------------------------"
            print *, "test_against_simple failed: Bder for case ", case
            print *, "Bder: ", bder
            print *, "SIMPLE: ", bder_ref(case, :)
            print *, "Relative error: ", abs(bder - bder_ref(case, :)) &
                /abs(bder_ref(case, :))
            print *, "Absolute error: ", abs(bder - bder_ref(case, :))
            test_failed = .true.
        end if
        hcovar = hcovar/cm2m
     if (not_same(hcovar, hcovar_ref(case, :), reltol_in=reltol, abstol_in=abstol)) then
            print *, "---------------------------------------------------"
            print *, "test_against_simple failed: hcovar for case ", case
            print *, "hcovar: ", hcovar
            print *, "SIMPLE: ", hcovar_ref(case, :)
            print *, "Relative error: ", abs(hcovar - hcovar_ref(case, :)) &
                /abs(hcovar_ref(case, :))
            print *, "Absolute error: ", abs(hcovar - hcovar_ref(case, :))
            test_failed = .true.
        end if
        hctrvr = hctrvr*cm2m
     if (not_same(hctrvr, hctrvr_ref(case, :), reltol_in=reltol, abstol_in=abstol)) then
            print *, "---------------------------------------------------"
            print *, "test_against_simple failed: hctrvr for case ", case
            print *, "hctrvr: ", hctrvr
            print *, "SIMPLE: ", hctrvr_ref(case, :)
            print *, "Relative error: ", abs(hctrvr - hctrvr_ref(case, :)) &
                /abs(hctrvr_ref(case, :))
            print *, "Absolute error: ", abs(hctrvr - hctrvr_ref(case, :))
            test_failed = .true.
        end if
        hcurl = hcurl*cm2m**2.0_dp
       if (not_same(hcurl, hcurl_ref(case, :), reltol_in=reltol, abstol_in=abstol)) then
            print *, "---------------------------------------------------"
            print *, "test_against_simple failed: hcurl for case ", case
            print *, "hcurl: ", hcurl
            print *, "SIMPLE: ", hcurl_ref(case, :)
            print *, "Relative error: ", abs(hcurl - hcurl_ref(case, :)) &
                /abs(hcurl_ref(case, :))
            print *, "Absolute error: ", abs(hcurl - hcurl_ref(case, :))
            test_failed = .true.
        end if
    end do

    if (test_failed) error stop

end program test_against_simple
