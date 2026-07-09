include(FetchContent)

# -DLIBNEO_PATH=<dir>: build against a local libneo checkout instead of fetching.
set(LIBNEO_PATH "" CACHE PATH "Local libneo source directory (leave empty to fetch)")

# -DLIBNEO_REF=<branch|tag|sha>: fetch this libneo revision instead of the pin.
set(LIBNEO_REF  "" CACHE STRING "libneo git ref (branch, tag, or sha) to fetch")

# rabe links only libneo's light targets (boozer, vmec_support, interpolate);
# skip its numerics-heavy core so BLAS/LAPACK and fortnum are not required here.
set(LIBNEO_BUILD_NUMERICS OFF CACHE BOOL "rabe needs only light libneo targets" FORCE)

# libneo is built in-tree from its own CMakeLists (either a local checkout
# or fetched), so it inherits our -Wall/-Wextra/... Debug flags. We do not
# care for internal warning in libneo as consumers so we append -w for its
# build only. As we fetch inside a function, this compiler flag change is
# constraint to the function scope and does not leak into the rest of the project.
function(rabe_add_libneo)
    string(APPEND CMAKE_Fortran_FLAGS_DEBUG   " -w")
    string(APPEND CMAKE_Fortran_FLAGS_RELEASE " -w")

    if(LIBNEO_PATH AND EXISTS "${LIBNEO_PATH}")
        message(STATUS "Using libneo in ${LIBNEO_PATH}")
        add_subdirectory("${LIBNEO_PATH}"
                         "${CMAKE_CURRENT_BINARY_DIR}/libneo"
                         EXCLUDE_FROM_ALL)
    else()
        if(LIBNEO_REF STREQUAL "")
            set(_libneo_ref "54e372f2ef15cf6cf0f36d734500c64d10cd5ada")
        else()
            set(_libneo_ref "${LIBNEO_REF}")
        endif()
        FetchContent_Declare(
            libneo
            GIT_REPOSITORY https://github.com/itpplasma/libneo.git
            GIT_TAG        ${_libneo_ref}
            EXCLUDE_FROM_ALL
        )
        FetchContent_MakeAvailable(libneo)
    endif()
endfunction()

rabe_add_libneo()
