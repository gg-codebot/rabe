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
