# End-to-end test of the electric_rabe restart file, see README.md
add_test(NAME RestartResume
         COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/run_restart_test.sh
                 ${CMAKE_BINARY_DIR}/electric_rabe
                 ${CMAKE_CURRENT_SOURCE_DIR}
                 ${CMAKE_CURRENT_SOURCE_DIR}/../integration/vmec/input/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(RestartResume PROPERTIES LABELS "slow")
