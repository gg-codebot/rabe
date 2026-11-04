# End-to-end test of the electric_rabe restart file: a run resumed from a
# restart file cut off mid-fieldline must reproduce an uninterrupted run exactly.
set(FIELD_FILE ${CMAKE_CURRENT_SOURCE_DIR}/../integration/vmec/input/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc)
set(FULL_DIR ${CMAKE_CURRENT_BINARY_DIR}/full)
set(RESUMED_DIR ${CMAKE_CURRENT_BINARY_DIR}/resumed)
foreach(run_dir ${FULL_DIR} ${RESUMED_DIR})
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/input/rabe.in ${FIELD_FILE}
         DESTINATION ${run_dir})
endforeach()

find_package(Python3 COMPONENTS Interpreter)
if(NOT Python3_FOUND)
    message(STATUS "Python3 not found - restart test disabled")
    return()
endif()

add_test(NAME RestartClean
         COMMAND ${CMAKE_COMMAND} -E rm -f
                 ${FULL_DIR}/electric_rabe.restart.dat
                 ${FULL_DIR}/electric_rabe.nc
                 ${RESUMED_DIR}/electric_rabe.restart.dat
                 ${RESUMED_DIR}/electric_rabe.nc)
set_tests_properties(RestartClean PROPERTIES FIXTURES_SETUP RestartCleaned)

add_test(NAME RestartFullRun
         COMMAND ${CMAKE_BINARY_DIR}/electric_rabe
         WORKING_DIRECTORY ${FULL_DIR})
set_tests_properties(RestartFullRun PROPERTIES
    FIXTURES_REQUIRED RestartCleaned
    FIXTURES_SETUP RestartFull)

add_test(NAME RestartTruncate
         COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/truncate_restart.py
                 ${FULL_DIR}/electric_rabe.restart.dat
                 ${RESUMED_DIR}/electric_rabe.restart.dat)
set_tests_properties(RestartTruncate PROPERTIES
    FIXTURES_REQUIRED RestartFull
    FIXTURES_SETUP RestartTruncated)

add_test(NAME RestartResumedRun
         COMMAND ${CMAKE_BINARY_DIR}/electric_rabe
         WORKING_DIRECTORY ${RESUMED_DIR})
set_tests_properties(RestartResumedRun PROPERTIES
    FIXTURES_REQUIRED RestartTruncated
    FIXTURES_SETUP RestartResumed
    PASS_REGULAR_EXPRESSION "Read +1 +of +2 +fieldlines")

add_test(NAME RestartCompareOutput
         COMMAND ${CMAKE_COMMAND} -E compare_files
                 ${FULL_DIR}/electric_rabe.nc
                 ${RESUMED_DIR}/electric_rabe.nc)
add_test(NAME RestartCompareRestartFile
         COMMAND ${CMAKE_COMMAND} -E compare_files
                 ${FULL_DIR}/electric_rabe.restart.dat
                 ${RESUMED_DIR}/electric_rabe.restart.dat)
set_tests_properties(RestartCompareOutput RestartCompareRestartFile PROPERTIES
    LABELS "slow"
    FIXTURES_REQUIRED RestartResumed)
