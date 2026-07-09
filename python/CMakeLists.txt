if(NOT SKBUILD_PROJECT_NAME)
    set(SKBUILD_PROJECT_NAME "rabe")
endif()

find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module NumPy)
find_program(F90WRAP_EXECUTABLE f90wrap REQUIRED)

set(F90WRAP_SOURCES
    ${PROJECT_SOURCE_DIR}/src/field/field_base.f90
    ${PROJECT_SOURCE_DIR}/src/field/fourier_field.f90
    ${PROJECT_SOURCE_DIR}/src/fieldline/fieldline.f90
    ${PROJECT_SOURCE_DIR}/src/vmec/boozer_field.f90
    ${PROJECT_SOURCE_DIR}/src/fieldline/surface_average.f90
    ${PROJECT_SOURCE_DIR}/src/fieldline/make_fieldline.f90
    ${PROJECT_SOURCE_DIR}/src/offset/deviation.f90
    ${PROJECT_SOURCE_DIR}/src/offset/coefficients.f90
    ${PROJECT_SOURCE_DIR}/src/error_handling/error_handling.f90
    ${PROJECT_SOURCE_DIR}/src/shaing_callen/shaing_callen.f90
)

foreach(src ${F90WRAP_SOURCES})
    get_filename_component(name ${src} NAME_WE)
    list(APPEND F90WRAP_GENERATED_F90
        ${CMAKE_CURRENT_BINARY_DIR}/f90wrap_${name}.f90)
endforeach()


add_custom_command(
    OUTPUT
        ${F90WRAP_GENERATED_F90}
        ${CMAKE_CURRENT_BINARY_DIR}/_rabe.c
    COMMAND ${F90WRAP_EXECUTABLE} --direct-c -P -M -C make_flock_of_fieldlines fourier_field_init -m rabe
            --f90-mod-name rabe._rabe
            -k ${CMAKE_CURRENT_SOURCE_DIR}/kind_map.py
            --joint-modules ${CMAKE_CURRENT_SOURCE_DIR}/joint_modules.py
            --class-names ${CMAKE_CURRENT_SOURCE_DIR}/class_names.py
            ${F90WRAP_SOURCES}
            --only
                fourier_field_t
                fourier_field_init
                boozer_field_t
                init_from_vmec
                fix_to_surface
                get_iota
                get_covariant_components
                field_t
                flock_of_fieldlines_t
                fieldline_t
                surface_average_t
                calc_surface_averages
                make_flock_of_fieldlines
                calc_deviation
                calc_nu_star_crit
                calc_finite_boundary_layer_correction
                calc_gradient_scaling_factor_r_eff
                calc_offset_coefficients
                set_unsafe_mode
                reset_failed_check_counter
                did_fail_any_sanity_check
                calc_lambda_lc
                get_non_omnigenous_remainder
    DEPENDS ${F90WRAP_SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/kind_map.py
            ${CMAKE_CURRENT_SOURCE_DIR}/joint_modules.py
            ${CMAKE_CURRENT_SOURCE_DIR}/class_names.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating f90wrap Python bindings"
    VERBATIM
)

Python3_add_library(_rabe MODULE WITH_SOABI
    ${F90WRAP_GENERATED_F90}
    ${CMAKE_CURRENT_BINARY_DIR}/_rabe.c
)

target_include_directories(_rabe PRIVATE ${Python3_NumPy_INCLUDE_DIRS})

target_link_libraries(_rabe PRIVATE
    rabe_lib
    offset_lib
    vmec_lib
    field_lib
    fieldline_lib
    utils_lib
    error_handling_lib
    shaing_callen_lib
    netcdf_lib
)

set_target_properties(_rabe PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

install(TARGETS _rabe DESTINATION "${SKBUILD_PROJECT_NAME}")
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rabe/
        DESTINATION "${SKBUILD_PROJECT_NAME}"
        FILES_MATCHING PATTERN "*.py")
