# Override with -DITP_FORTIO_REF=<branch|tag|sha>; see README.
set(RABE_FORTIO_REF 11afd0bd1af0c99ea4e9a1c0df683dbdcc299b69)
itp_add_dependency(fortio REF ${RABE_FORTIO_REF}
    GIT_REPOSITORY https://github.com/lazy-fortran/fortio.git
    TARGET fortio)
