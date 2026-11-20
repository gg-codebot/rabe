# Override with -DITP_LIBNEO_REF=<branch|tag|sha> or use a local checkout with
# -DFETCHCONTENT_SOURCE_DIR_LIBNEO=<dir>; see README.
set(RABE_LIBNEO_REF 7262c32bf0406594a8ef27e84654d249aef2cb86)

# libneo is built in-tree from its own CMakeLists (either a local checkout
# or fetched), so it inherits our -Wall/-Wextra/... Debug flags. We do not
# care for internal warning in libneo as consumers so we append -w for its
# build only. As we fetch inside a function, this compiler flag change is
# constraint to the function scope and does not leak into the rest of the project.
function(rabe_add_libneo)
    string(APPEND CMAKE_Fortran_FLAGS_DEBUG   " -w")
    string(APPEND CMAKE_Fortran_FLAGS_RELEASE " -w")

    # rabe links only libneo's light targets (boozer, vmec_support,
    # interpolate); skip its numerics-heavy core so BLAS/LAPACK and fortnum
    # are not required here.
    itp_add_dependency(libneo REF ${RABE_LIBNEO_REF} TARGET LIBNEO::boozer
        OPTIONS "LIBNEO_BUILD_NUMERICS OFF")
endfunction()

rabe_add_libneo()
