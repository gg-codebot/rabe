add_library(netcdf_lib STATIC
    netcdf.f90
  )

target_link_libraries(netcdf_lib PRIVATE utils_lib fortio::fortio)
target_include_directories(netcdf_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
