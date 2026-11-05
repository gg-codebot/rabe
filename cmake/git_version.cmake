execute_process(
    COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${SOURCE_DIR}
    OUTPUT_VARIABLE GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)

if(NOT GIT_HASH)
    set(GIT_HASH "unknown")
else()
    execute_process(
        COMMAND git status --porcelain --untracked-files=no
        WORKING_DIRECTORY ${SOURCE_DIR}
        OUTPUT_VARIABLE GIT_CHANGES
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(GIT_CHANGES)
        set(GIT_HASH "${GIT_HASH}-dirty")
    endif()
endif()

configure_file(
    ${SOURCE_DIR}/src/git_version.f90.in
    ${OUTPUT_FILE}
    @ONLY
)
