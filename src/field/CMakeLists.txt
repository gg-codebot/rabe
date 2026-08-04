add_library(field_lib STATIC field_base.f90)
target_link_libraries(field_lib PRIVATE utils_lib)
target_include_directories(field_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

add_library(fourier_field_lib STATIC fourier_field.f90)
# interpolate PUBLIC: fourier_field_t's public component `spl` is
# type(SplineData2D) from interpolate, so consumers declaring/using
# fourier_field_t need interpolate.mod transitively resolvable too.
target_link_libraries(fourier_field_lib PRIVATE field_lib utils_lib PUBLIC interpolate)
target_include_directories(fourier_field_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
