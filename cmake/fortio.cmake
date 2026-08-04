include(FetchContent)

if(NOT TARGET fortio)
    FetchContent_Declare(
        fortio
        GIT_REPOSITORY https://github.com/lazy-fortran/fortio.git
        GIT_TAG ce5c7257563648dc5afbdfbe5b5a91181bf06912
    )
    FetchContent_MakeAvailable(fortio)
endif()
