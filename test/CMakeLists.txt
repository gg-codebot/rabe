include(CheckLinkerFlag)

# A few test programs pass internal procedures as callback arguments (e.g. to
# quadpack via integrate_1d), which requires an executable-stack trampoline.
# this is not a safety issue for test executables, and infrastructures that
# do not allow trampolines a priori will fail the build with an error anyway.
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    # gfortran already marks the stack executable (GNU_STACK RWE) for this
    # pattern; we only need to silence the resulting linker warning.
    check_linker_flag(Fortran "-Wl,--no-warn-execstack" SUPPORTS_NO_WARN_EXECSTACK)
    if(SUPPORTS_NO_WARN_EXECSTACK)
        add_link_options(-Wl,--no-warn-execstack)
    endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "LLVMFlang")
    # Unlike gfortran, flang does NOT mark the stack executable for this same
    # trampoline pattern (GNU_STACK comes out RW), so request it explicitly.
    check_linker_flag(Fortran "-Wl,-z,execstack" SUPPORTS_EXECSTACK)
    if(SUPPORTS_EXECSTACK)
        add_link_options(-Wl,-z,execstack)
    endif()
endif()

add_subdirectory(plot_lib)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/plot_lib)

# Under flang, CMake/Ninja's dyndep dependency
# scanning always writes a preprocessed copy to disk before compiling when
# preprocessing is enabled, and flang's preprocessor output corrupts
# backslash escapes in string literals (e.g. LaTeX labels like "$\pi$"
# become "$\\pi$"), breaking matplotlib label parsing at runtime. Since
# nothing here needs preprocessing anyway, skip it outright to avoid the bug.
set(CMAKE_Fortran_PREPROCESS OFF)

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
