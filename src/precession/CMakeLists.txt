add_library(precession_lib STATIC
  grid.f90
  splines_instance.f90
  precession.f90
  bounce.f90
  orbit.f90
  params.f90
  field_instance.f90
  odeint_allroutines.f
  read_precession.f90
  )

include(libneo)
target_link_libraries(precession_lib PRIVATE utils_lib
                      field_lib
                      fieldline_lib
                      vmec_lib
                      LIBNEO::interpolate
                      LIBNEO::odeint
                      LIBNEO::util
                      LIBNEO::vmec_support)
target_include_directories(precession_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
