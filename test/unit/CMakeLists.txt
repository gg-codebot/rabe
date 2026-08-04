add_executable(test_handle_multiple_maxima.x test_handle_multiple_maxima.f90)
target_link_libraries(test_handle_multiple_maxima.x utils_lib fieldline_lib)
add_test(TestHandleMultipleMaxima test_handle_multiple_maxima.x)
set_tests_properties(TestHandleMultipleMaxima PROPERTIES LABELS "quick")

add_executable(test_find_extrema.x test_find_extrema.f90)
target_link_libraries(test_find_extrema.x utils_lib)
add_test(TestFindExtrema test_find_extrema.x)
set_tests_properties(TestFindExtrema PROPERTIES LABELS "quick")

add_executable(test_integrate_1d.x test_integrate_1d.f90)
target_link_libraries(test_integrate_1d.x utils_lib)
add_test(TestIntegrate1d test_integrate_1d.x)
set_tests_properties(TestIntegrate1d PROPERTIES LABELS "quick")

add_executable(test_integrate_1d_substituted.x test_integrate_1d_substituted.f90)
target_link_libraries(test_integrate_1d_substituted.x utils_lib)
add_test(TestIntegrate1dSubstituted test_integrate_1d_substituted.x)
set_tests_properties(TestIntegrate1dSubstituted PROPERTIES LABELS "quick")

add_executable(test_sum_trapez_1d.x test_sum_trapez_1d.f90)
target_link_libraries(test_sum_trapez_1d.x utils_lib)
add_test(TestISumTrapez test_sum_trapez_1d.x)
set_tests_properties(TestISumTrapez PROPERTIES LABELS "quick")

add_executable(test_integrate_wrapper.x test_integrate_wrapper.f90)
target_link_libraries(test_integrate_wrapper.x utils_lib)
add_test(TestIntegrateWrapper test_integrate_wrapper.x)
set_tests_properties(TestIntegrateWrapper PROPERTIES LABELS "quick")

add_executable(test_real_ft.x test_real_ft.f90)
target_link_libraries(test_real_ft.x utils_lib)
add_test(TestRealFt test_real_ft.x)
set_tests_properties(TestRealFt PROPERTIES LABELS "quick")

add_executable(test_fit_functions.x test_fit_functions.f90)
target_link_libraries(test_fit_functions.x offset_lib utils_lib)
add_test(TestFitFunctions test_fit_functions.x)
set_tests_properties(TestFitFunctions PROPERTIES LABELS "quick")

add_executable(test_netcdf.x test_netcdf.f90)
target_link_libraries(test_netcdf.x utils_lib netcdf_lib)
add_test(TestNetcdf test_netcdf.x)
set_tests_properties(TestNetcdf PROPERTIES LABELS "quick")

add_executable(test_integrate_over_eta_grid.x test_integrate_over_eta_grid.f90)
target_link_libraries(test_integrate_over_eta_grid.x utils_lib shaing_callen_lib)
add_test(TestIntegrateOverEtaGrid test_integrate_over_eta_grid.x)
set_tests_properties(TestIntegrateOverEtaGrid PROPERTIES LABELS "quick")

add_executable(test_integrate_over_phi_grid.x test_integrate_over_phi_grid.f90)
target_link_libraries(test_integrate_over_phi_grid.x utils_lib shaing_callen_lib fieldline_lib)
add_test(TestIntegrateOverPhiGrid test_integrate_over_phi_grid.x)
set_tests_properties(TestIntegrateOverPhiGrid PROPERTIES LABELS "quick")

add_executable(test_calc_periodic_dydx.x test_calc_periodic_dydx.f90)
target_link_libraries(test_calc_periodic_dydx.x utils_lib shaing_callen_lib)
add_test(TestCalcPeriodicDerivative test_calc_periodic_dydx.x)
set_tests_properties(TestCalcPeriodicDerivative PROPERTIES LABELS "quick")

add_executable(test_read_namelist.x test_read_namelist.f90)
target_link_libraries(test_read_namelist.x read_file_lib utils_lib)
add_test(TestReadNamelist test_read_namelist.x)
set_tests_properties(TestReadNamelist PROPERTIES LABELS "quick")

add_executable(test_linspace.x test_linspace.f90)
target_link_libraries(test_linspace.x utils_lib)
add_test(TestLinspace test_linspace.x)
set_tests_properties(TestLinspace PROPERTIES LABELS "quick")

add_executable(test_diophantine.x test_diophantine.f90)
target_link_libraries(test_diophantine.x utils_lib)
add_test(TestDiophantine test_diophantine.x)
set_tests_properties(TestDiophantine PROPERTIES LABELS "quick")

add_executable(test_error_handling.x test_error_handling.f90)
target_link_libraries(test_error_handling.x error_handling_lib)
add_test(TestErrorHandling test_error_handling.x)
set_tests_properties(TestErrorHandling PROPERTIES LABELS "quick")

add_executable(plot_fit_functions.x plot_fit_functions.f90)
target_link_libraries(plot_fit_functions.x offset_lib utils_lib myplot_lib)
add_test(PlotFitFunctions plot_fit_functions.x)
set_tests_properties(PlotFitFunctions PROPERTIES LABELS "plot")
