file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/input DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../vmec/input DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_trace_till_bounce.x test_trace_till_bounce.f90)
target_link_libraries(test_trace_till_bounce.x rabe_lib helpers_lib)
add_test(TestTraceTillBounce test_trace_till_bounce.x)
set_tests_properties(TestTraceTillBounce PROPERTIES LABELS "quick")

add_executable(test_precession_analytic.x test_precession_analytic.f90)
target_link_libraries(test_precession_analytic.x rabe_lib helpers_lib myplot_lib)
add_test(TestPrecessionAnalytic test_precession_analytic.x)
set_tests_properties(TestPrecessionAnalytic PROPERTIES LABELS "current")

add_executable(test_bounce_integrals.x test_bounce_integrals.f90)
target_link_libraries(test_bounce_integrals.x rabe_lib helpers_lib myplot_lib)
add_test(TestBounceIntegrals test_bounce_integrals.x)
set_tests_properties(TestBounceIntegrals PROPERTIES LABELS "slow")

add_executable(test_radial_drift.x test_radial_drift.f90)
target_link_libraries(test_radial_drift.x rabe_lib helpers_lib myplot_lib)
add_test(TestRadialDrift test_radial_drift.x)
set_tests_properties(TestRadialDrift PROPERTIES LABELS "slow")

add_executable(test_zero_precession.x test_zero_precession.f90)
target_link_libraries(test_zero_precession.x rabe_lib helpers_lib)
add_test(TestZeroPrecession test_zero_precession.x)
set_tests_properties(TestZeroPrecession PROPERTIES LABELS "slow")

add_executable(test_g_zero_precession.x test_g_zero_precession.f90)
target_link_libraries(test_g_zero_precession.x rabe_lib helpers_lib)
add_test(TestGZeroPrecession test_g_zero_precession.x)
set_tests_properties(TestGZeroPrecession PROPERTIES LABELS "slow")
