add_executable(test_plot.x test_plot.f90)
target_link_libraries(test_plot.x myplot_lib utils_lib)
add_test(TestPlot test_plot.x)
set_tests_properties(TestPlot PROPERTIES LABELS "plot")

add_executable(plot_anti_sigma.x plot_anti_sigma.f90)
target_link_libraries(plot_anti_sigma.x myplot_lib
                                             helpers_lib
                                             utils_lib
                                             fieldline_lib)
add_test(PlotAntiSigma plot_anti_sigma.x)
set_tests_properties(PlotAntiSigma PROPERTIES LABELS "plot")

add_executable(plot_deviation_helical_anti_sigma.x plot_deviation_helical_anti_sigma.f90)
target_link_libraries(plot_deviation_helical_anti_sigma.x myplot_lib
                                                          offset_lib
                                                          helpers_lib
                                                          utils_lib
                                                          fieldline_lib)
add_test(PlotDeviationHelicalAntiSigma plot_deviation_helical_anti_sigma.x)
set_tests_properties(PlotDeviationHelicalAntiSigma PROPERTIES LABELS "plot")

add_executable(plot_helical_anti_sigma_small_aspect.x plot_helical_anti_sigma_small_aspect.f90)
target_link_libraries(plot_helical_anti_sigma_small_aspect.x myplot_lib
                                                             offset_lib
                                                             helpers_lib
                                                             utils_lib
                                                             fieldline_lib)
add_test(PlotHelicalAntiSigmaSmallAspect plot_helical_anti_sigma_small_aspect.x)
set_tests_properties(PlotHelicalAntiSigmaSmallAspect PROPERTIES LABELS "plot")

add_executable(plot_deviation_poloidal_anti_sigma.x plot_deviation_poloidal_anti_sigma.f90)
target_link_libraries(plot_deviation_poloidal_anti_sigma.x myplot_lib
                                                          offset_lib
                                                          helpers_lib
                                                          utils_lib
                                                          fieldline_lib)
add_test(PlotDeviationPoloidalAntiSigma plot_deviation_poloidal_anti_sigma.x)
set_tests_properties(PlotDeviationPoloidalAntiSigma PROPERTIES LABELS "plot")

add_executable(plot_local_drift.x plot_local_drift.f90)
target_link_libraries(plot_local_drift.x myplot_lib
                                             helpers_lib
                                             utils_lib
                                             fieldline_lib)
add_test(PlotLocalDrift plot_local_drift.x)
set_tests_properties(PlotLocalDrift PROPERTIES LABELS "plot")

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../integration/vmec/input
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(plot_fourier_field_timing.x plot_fourier_field_timing.f90)
target_link_libraries(plot_fourier_field_timing.x myplot_lib rabe_lib helpers_lib)
add_test(PlotFourierFieldTiming plot_fourier_field_timing.x)
set_tests_properties(PlotFourierFieldTiming PROPERTIES LABELS "plot")

add_executable(plot_fourier_field_accuracy.x plot_fourier_field_accuracy.f90)
target_link_libraries(plot_fourier_field_accuracy.x myplot_lib rabe_lib helpers_lib)
add_test(PlotFourierFieldAccuracy plot_fourier_field_accuracy.x)
set_tests_properties(PlotFourierFieldAccuracy PROPERTIES LABELS "plot")

add_executable(plot_bounce_time_averages.x plot_bounce_time_averages.f90)
target_link_libraries(plot_bounce_time_averages.x myplot_lib
                                                     rabe_lib
                                                     helpers_lib)
add_test(PlotBounceTimeAverages plot_bounce_time_averages.x)
set_tests_properties(PlotBounceTimeAverages PROPERTIES LABELS "plot")

add_executable(plot_bounce_orbits.x plot_bounce_orbits.f90)
target_link_libraries(plot_bounce_orbits.x myplot_lib
                                                 rabe_lib
                                                 helpers_lib)
add_test(PlotBounceOrbits plot_bounce_orbits.x)
set_tests_properties(PlotBounceOrbits PROPERTIES LABELS "plot")

add_executable(plot_integrate_radial_drift.x plot_integrate_radial_drift.f90)
target_link_libraries(plot_integrate_radial_drift.x rabe_lib
                                                   helpers_lib)
add_test(PlotIntegrateRadialDrift plot_integrate_radial_drift.x)
set_tests_properties(PlotIntegrateRadialDrift PROPERTIES LABELS "plot")
