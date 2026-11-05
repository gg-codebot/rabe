cmake_minimum_required(VERSION 3.24)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(rabe)
enable_language(Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
include(compiler)

# required for fetching dependencies
find_package(Git REQUIRED)

include(fortio)

# runs on every build to detect uncommitted changes
add_custom_target(git_version_target
    COMMAND ${CMAKE_COMMAND}
        -DSOURCE_DIR=${PROJECT_SOURCE_DIR}
        -DOUTPUT_FILE=${PROJECT_BINARY_DIR}/src/git_version.f90
        -P ${PROJECT_SOURCE_DIR}/cmake/git_version.cmake
    BYPRODUCTS ${PROJECT_BINARY_DIR}/src/git_version.f90
    COMMENT "Generating git_version.f90 with current git hash"
    VERBATIM
)

add_subdirectory(src) # for rabe_lib

add_executable(${PROJECT_NAME} app/main.f90)
target_link_libraries(${PROJECT_NAME} rabe_lib)
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "rabe.x")

add_executable(electric_rabe app/electric_rabe.f90)
target_link_libraries(electric_rabe rabe_lib)

enable_testing()
add_subdirectory(test)

option(BUILD_PYTHON_BINDINGS "Build Python bindings via f90wrap" OFF)
if(BUILD_PYTHON_BINDINGS)
    add_subdirectory(python)
endif()
