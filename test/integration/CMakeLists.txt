add_executable(test_fieldline.x test_fieldline.f90)
target_link_libraries(test_fieldline.x helpers_lib utils_lib fieldline_lib)
add_test(TestFieldline test_fieldline.x)
set_tests_properties(TestFieldline PROPERTIES LABELS "quick")

add_executable(test_make_flock_of_fieldlines.x test_make_flock_of_fieldlines.f90)
target_link_libraries(test_make_flock_of_fieldlines.x helpers_lib utils_lib fieldline_lib)
add_test(TestMakeFlockOfFieldlines test_make_flock_of_fieldlines.x)
set_tests_properties(TestMakeFlockOfFieldlines PROPERTIES LABELS "quick")

add_executable(test_integration_along_fieldline.x test_integration_along_fieldline.f90)
target_link_libraries(test_integration_along_fieldline.x helpers_lib utils_lib fieldline_lib)
add_test(TestIntegrationFieldline test_integration_along_fieldline.x)
set_tests_properties(TestIntegrationFieldline PROPERTIES LABELS "quick")

add_executable(test_find_global_maximum.x test_find_global_maximum.f90)
target_link_libraries(test_find_global_maximum.x helpers_lib utils_lib fieldline_lib)
add_test(TestGlobalMaximum test_find_global_maximum.x)
set_tests_properties(TestGlobalMaximum PROPERTIES LABELS "quick")

add_executable(test_find_turning_points.x test_find_turning_points.f90)
target_link_libraries(test_find_turning_points.x rabe_lib helpers_lib)
add_test(TestFindTurningPoints test_find_turning_points.x)
set_tests_properties(TestFindTurningPoints PROPERTIES LABELS "quick")

add_executable(test_anti_sigma_analytic.x test_anti_sigma_analytic.f90)
target_link_libraries(test_anti_sigma_analytic.x helpers_lib utils_lib fieldline_lib)
add_test(TestAntiSigmaAnalytic test_anti_sigma_analytic.x)
set_tests_properties(TestAntiSigmaAnalytic PROPERTIES LABELS "quick")

add_executable(test_pert_anti_sigma_analytic.x test_pert_anti_sigma_analytic.f90)
target_link_libraries(test_pert_anti_sigma_analytic.x helpers_lib utils_lib fieldline_lib)
add_test(TestPertAntiSigmaAnalytic test_pert_anti_sigma_analytic.x)
set_tests_properties(TestPertAntiSigmaAnalytic PROPERTIES LABELS "quick")

add_executable(test_I_integral_helical_anti_sigma.x test_I_integral_helical_anti_sigma.f90)
target_link_libraries(test_I_integral_helical_anti_sigma.x offset_lib
                                                           helpers_lib
                                                           myplot_lib
                                                           utils_lib
                                                           fieldline_lib)
add_test(TestIIntegralHelicalAntiSigma test_I_integral_helical_anti_sigma.x)
set_tests_properties(TestIIntegralHelicalAntiSigma PROPERTIES LABELS "quick")

add_executable(test_lambda_helical_anti_sigma.x test_lambda_helical_anti_sigma.f90)
target_link_libraries(test_lambda_helical_anti_sigma.x helpers_lib
                                                       myplot_lib
                                                       utils_lib
                                                       fieldline_lib)
add_test(TestLambdaHelicalAntiSigma test_lambda_helical_anti_sigma.x)
set_tests_properties(TestLambdaHelicalAntiSigma PROPERTIES LABELS "quick")


add_executable(test_set_fieldline_labels.x test_set_fieldline_labels.f90)
target_link_libraries(test_set_fieldline_labels.x offset_lib helpers_lib utils_lib fieldline_lib)
add_test(TestSetFieldlineLabels test_set_fieldline_labels.x)
set_tests_properties(TestSetFieldlineLabels PROPERTIES LABELS "quick")

add_executable(test_check_field_origin.x test_check_field_origin.f90)
target_link_libraries(test_check_field_origin.x helpers_lib utils_lib fieldline_lib)
add_test(TestCheckFieldOrigin test_check_field_origin.x)
set_tests_properties(TestCheckFieldOrigin PROPERTIES LABELS "quick")

add_executable(test_shaing_callen_against_quasi_symmetric.x test_shaing_callen_against_quasi_symmetric.f90)
target_link_libraries(test_shaing_callen_against_quasi_symmetric.x helpers_lib myplot_lib
                                                                    utils_lib fieldline_lib shaing_callen_lib)
add_test(TestShaingCallenQS test_shaing_callen_against_quasi_symmetric.x)
set_tests_properties(TestShaingCallenQS PROPERTIES LABELS "quick")

add_subdirectory(vmec)
add_subdirectory(precession)
