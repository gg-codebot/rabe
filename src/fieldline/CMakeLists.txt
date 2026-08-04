add_library(fieldline_lib STATIC
  fieldline.f90
  fieldline_labels.f90
  field_along_fieldline.f90
  field_checks.f90
  make_fieldline.f90
  fieldline_integrals.f90
  fieldline_integrands.f90
  surface_average.f90
  )

# field_lib PUBLIC: flang (unlike gfortran) requires field_base.mod to be
# transitively resolvable wherever a type/procedure from it is used, not just
# within fieldline_lib itself (e.g. make_fieldline.f90 exposes field_t).
target_link_libraries(fieldline_lib PRIVATE utils_lib error_handling_lib
  PUBLIC field_lib)
target_include_directories(fieldline_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
