# electric_rabe reads rabe.in from its working directory, so it runs in its own
# build directory
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/input/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${PROJECT_SOURCE_DIR}/test/integration/vmec/input/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# start without saved bounce integrals from an earlier run
add_test(NAME GoldenElectricRun
         COMMAND sh -c "rm -f saved_bounce_integrals.dat && ${CMAKE_BINARY_DIR}/electric_rabe"
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(GoldenElectricRun PROPERTIES
    FIXTURES_SETUP GoldenElectricOutput)

if(Python3_FOUND)
    add_test(NAME GoldenElectricCompare
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../compare_golden.py
                     ${CMAKE_CURRENT_SOURCE_DIR}/expected/electric_rabe.nc
                     ${CMAKE_CURRENT_BINARY_DIR}/electric_rabe.nc)
    set_tests_properties(GoldenElectricCompare PROPERTIES
        LABELS "golden"
        FIXTURES_REQUIRED GoldenElectricOutput
        SKIP_RETURN_CODE 77)
endif()
