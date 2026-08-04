import rabe

# Top-level submodules
for mod in [
    "boozer_field",
    "fourier_field",
    "fieldline_mod",
    "surface_average_mod",
    "error_handling",
    "coefficients",
]:
    assert hasattr(rabe, mod), f"Missing submodule: {mod}"

# Classes
assert hasattr(
    rabe.boozer_field, "BoozerField"
), "Missing rabe.boozer_field.BoozerField"
assert hasattr(
    rabe.fourier_field, "FourierField"
), "Missing rabe.fourier_field.FourierField"
assert hasattr(
    rabe.fieldline_mod, "FlockOfFieldlines"
), "Missing rabe.fieldline_mod.FlockOfFieldlines"
assert hasattr(
    rabe.fieldline_mod, "fieldline_t"
), "Missing rabe.fieldline_mod.fieldline_t"
assert hasattr(
    rabe.surface_average_mod, "SurfaceAverage"
), "Missing rabe.surface_average_mod.SurfaceAverage"

# Methods on FlockOfFieldlines
flock_methods = [
    "calc_surface_averages",
    "calc_deviation",
    "calc_nu_star_crit",
    "calc_finite_boundary_layer_correction",
    "calc_gradient_scaling_factor_r_eff",
    "calc_offset_coefficients",
    "calc_lambda_lc",
    "get_non_omnigenous_remainder",
]
for m in flock_methods:
    assert hasattr(
        rabe.fieldline_mod.FlockOfFieldlines, m
    ), f"Missing rabe.fieldline_mod.FlockOfFieldlines.{m}"

# Methods on BoozerField
boozer_methods = [
    "init_from_vmec",
    "fix_to_surface",
    "get_iota",
]
for m in boozer_methods:
    assert hasattr(
        rabe.boozer_field.BoozerField, m
    ), f"Missing rabe.boozer_field.BoozerField.{m}"

# Standalone functions
error_handling_fns = [
    "set_unsafe_mode",
    "reset_failed_check_counter",
    "did_fail_any_sanity_check",
]
for fn in error_handling_fns:
    assert hasattr(rabe.error_handling, fn), f"Missing rabe.error_handling.{fn}"

print("All symbols present.")
