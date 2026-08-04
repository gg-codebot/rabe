add_library(myplot_lib STATIC
  myplot_module.f90
  )

add_subdirectory(pyplot-fortran)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/pyplot-fortran)
# pyplot-fortran (vendored) doesn't declare its own .mod include dir as a
# target property, only via the directory-scoped include_directories() above,
# which doesn't propagate to consumers in other directories. Declare it here
# so flang's stricter transitive-.mod resolution (see PUBLIC link below) can
# actually find pyplot_module.mod from outside this directory.
target_include_directories(pyplot_fortran PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/pyplot-fortran)
# PUBLIC: flang (unlike gfortran) requires pyplot_module.mod to be transitively
# resolvable wherever myplot_module is used, since myplot_module's public
# interface (myplot_module.mod) re-exposes the pyplot derived type.
target_link_libraries(myplot_lib PUBLIC pyplot_fortran)
target_link_libraries(myplot_lib PRIVATE utils_lib)
