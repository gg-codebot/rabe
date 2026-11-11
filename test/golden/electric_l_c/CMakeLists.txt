# same case as ../electric, but the collisionality grid is given as l_c instead
# of nu_star; the result must match the electric golden record
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/input/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${PROJECT_SOURCE_DIR}/test/integration/vmec/input/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME GoldenElectricLcRun
         COMMAND sh -c "rm -f saved_bounce_integrals.dat && ${CMAKE_BINARY_DIR}/electric_rabe"
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(GoldenElectricLcRun PROPERTIES
    FIXTURES_SETUP GoldenElectricLcOutput)

if(Python3_FOUND)
    add_test(NAME GoldenElectricLcCompare
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../compare_golden.py
                     ${CMAKE_CURRENT_SOURCE_DIR}/../electric/expected/electric_rabe.nc
                     ${CMAKE_CURRENT_BINARY_DIR}/electric_rabe.nc)
    set_tests_properties(GoldenElectricLcCompare PROPERTIES
        LABELS "golden"
        FIXTURES_REQUIRED GoldenElectricLcOutput
        SKIP_RETURN_CODE 77)
endif()
