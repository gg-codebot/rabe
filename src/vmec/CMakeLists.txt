add_library(vmec_lib STATIC
            boozer_field.f90
            )

include(libneo)
target_link_libraries(vmec_lib PUBLIC utils_lib field_lib)
target_link_libraries(vmec_lib PRIVATE
                      LIBNEO::boozer
                      LIBNEO::vmec_support)
target_include_directories(vmec_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
