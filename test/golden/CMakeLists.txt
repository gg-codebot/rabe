file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/input/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../integration/vmec/input/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/expected DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/compare_golden.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME GoldenRecordRun
         COMMAND ${CMAKE_BINARY_DIR}/rabe.x
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(GoldenRecordRun PROPERTIES
    FIXTURES_SETUP GoldenOutput)

find_package(Python3 COMPONENTS Interpreter)
if(Python3_FOUND)
    add_test(NAME GoldenRecordCompare
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/compare_golden.py
                     ${CMAKE_CURRENT_BINARY_DIR}/expected/rabe.nc
                     ${CMAKE_CURRENT_BINARY_DIR}/rabe.nc
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(GoldenRecordCompare PROPERTIES
        LABELS "golden"
        FIXTURES_REQUIRED GoldenOutput
        SKIP_RETURN_CODE 77)
else()
    message(STATUS "Python3 not found - golden record comparison test disabled")
endif()

# electric_rabe: separate directory, it reads its own rabe.in
set(ELECTRIC_DIR ${CMAKE_CURRENT_BINARY_DIR}/electric)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/electric_input/ DESTINATION ${ELECTRIC_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../integration/vmec/input/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc DESTINATION ${ELECTRIC_DIR})

# start without saved bounce integrals from an earlier run
add_test(NAME GoldenElectricRun
         COMMAND sh -c "rm -f saved_bounce_integrals.dat && ${CMAKE_BINARY_DIR}/electric_rabe"
         WORKING_DIRECTORY ${ELECTRIC_DIR})
set_tests_properties(GoldenElectricRun PROPERTIES
    FIXTURES_SETUP GoldenElectricOutput)

if(Python3_FOUND)
    add_test(NAME GoldenElectricCompare
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/compare_golden.py
                     ${CMAKE_CURRENT_BINARY_DIR}/expected/electric_rabe.nc
                     ${ELECTRIC_DIR}/electric_rabe.nc)
    set_tests_properties(GoldenElectricCompare PROPERTIES
        LABELS "golden"
        FIXTURES_REQUIRED GoldenElectricOutput
        SKIP_RETURN_CODE 77)
endif()
