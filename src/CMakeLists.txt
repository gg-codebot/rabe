add_library(git_version_lib STATIC
  ${CMAKE_CURRENT_BINARY_DIR}/git_version.f90
)

add_dependencies(git_version_lib git_version_target)
target_include_directories(git_version_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

add_library(read_file_lib STATIC read_file.f90)
target_include_directories(read_file_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(read_file_lib PRIVATE utils_lib fieldline_lib)

add_subdirectory(error_handling)
add_subdirectory(utils)
add_subdirectory(netcdf)
add_subdirectory(vmec) # also provides the `interpolate` target via include(libneo)
add_subdirectory(field)
add_subdirectory(fieldline)
add_subdirectory(offset)
add_subdirectory(shaing_callen)

add_library(rabe_lib INTERFACE)
target_link_libraries(rabe_lib INTERFACE
    git_version_lib
    read_file_lib
    utils_lib
    error_handling_lib
    field_lib
    fourier_field_lib
    fieldline_lib
    shaing_callen_lib
    offset_lib
    netcdf_lib
    vmec_lib
)
