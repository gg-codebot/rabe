add_library(shaing_callen_lib STATIC
  shaing_callen.f90
  shaing_callen_integration.f90
  shaing_callen_wrappers.f90
  )

# fieldline_lib PUBLIC: flang (unlike gfortran) requires fieldline_mod.mod to be
# transitively resolvable wherever a type/procedure from it is used, not just
# within shaing_callen_lib itself (e.g. shaing_callen_integration.f90 uses
# fieldline_t from fieldline_mod).
target_link_libraries(shaing_callen_lib PRIVATE utils_lib field_lib
  PUBLIC fieldline_lib)
target_include_directories(shaing_callen_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
