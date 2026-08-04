cmake_minimum_required(VERSION 3.24)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(rabe)
enable_language(Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
include(compiler)

# required for fetching dependencies
find_package(Git REQUIRED)

include(fortio)

if(EXISTS ${PROJECT_SOURCE_DIR}/.git/HEAD)
    set(GIT_HEAD_DEPS
        ${PROJECT_SOURCE_DIR}/.git/HEAD
        ${PROJECT_SOURCE_DIR}/.git/index
    )
endif()

add_custom_command(
    OUTPUT ${PROJECT_BINARY_DIR}/src/git_version.f90
    COMMAND ${CMAKE_COMMAND}
        -DSOURCE_DIR=${PROJECT_SOURCE_DIR}
        -DOUTPUT_FILE=${PROJECT_BINARY_DIR}/src/git_version.f90
        -P ${PROJECT_SOURCE_DIR}/cmake/git_version.cmake
    DEPENDS
        ${PROJECT_SOURCE_DIR}/src/git_version.f90.in
        ${GIT_HEAD_DEPS}
    COMMENT "Generating git_version.f90 with current git hash"
    VERBATIM
)

add_custom_target(git_version_target
    DEPENDS ${PROJECT_BINARY_DIR}/src/git_version.f90
)

add_subdirectory(src) # for rabe_lib

add_executable(${PROJECT_NAME} app/main.f90)
target_link_libraries(${PROJECT_NAME} rabe_lib)
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "rabe.x")

enable_testing()
add_subdirectory(test)

option(BUILD_PYTHON_BINDINGS "Build Python bindings via f90wrap" OFF)
if(BUILD_PYTHON_BINDINGS)
    add_subdirectory(python)
endif()
