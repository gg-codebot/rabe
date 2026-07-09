if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    # A few test programs pass internal procedures as callback arguments,
    # which gfortran implements via stack trampolines.
    # We silence respective warnings in test/ only as for the executables
    # this is not a safety issue and infrastructures that do not allow
    # trampolins a priori will fail the build with an error anyway.
    # --no-warn-execstack is a fairly recent GNU ld addition, so probe first
    include(CheckLinkerFlag)
    check_linker_flag(Fortran "-Wl,--no-warn-execstack" SUPPORTS_NO_WARN_EXECSTACK)
    if(SUPPORTS_NO_WARN_EXECSTACK)
        add_link_options(-Wl,--no-warn-execstack)
    endif()
endif()

add_subdirectory(plot_lib)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/plot_lib)

add_subdirectory(helpers)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/helpers)

enable_testing()
add_subdirectory(unit)
add_subdirectory(integration)

find_package(Python3 COMPONENTS Interpreter)
if(Python3_FOUND)
    add_subdirectory(plot)
else()
    message(STATUS "Python3 not found - plot tests disabled")
endif()
add_subdirectory(golden)

 # Private test data and comparisons (not available in public repo)
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/internal/CMakeLists.txt)
    add_subdirectory(internal)
endif()
