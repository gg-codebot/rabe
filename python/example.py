from pathlib import Path

from rabe.boozer_field import BoozerField
from rabe.fieldline_mod import FlockOfFieldlines

FIELD_FILE = str(
    Path(__file__).parent
    / "../test/integration/vmec/input/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc"
)
print(f"Using field file: {FIELD_FILE}")
M_POL = -1.0
N_TOR = 4.0
S_TOR = [0.255102041, 0.316326531, 0.5, 0.704081633, 0.908163265]
SIGN_SQRTG = -1
MAX_N_FIELDLINES = 200
SHOULD_CALC_SHAING_CALLEN = True
N_ETA = 100

field = BoozerField()
field.init_from_vmec(
    FIELD_FILE, radial_spline_order=5, angular_spline_order=5, grid_refinement=6
)
print(f"R = {field.r}")

for s in S_TOR:
    field.fix_to_surface(s)
    iota = field.get_iota(s)
    nfp = field.nfp

    flock = FlockOfFieldlines(MAX_N_FIELDLINES, iota, field, M_POL, N_TOR, nfp)

    dr_datheta = flock.calc_gradient_scaling_factor_r_eff(
        field.psi_tor_edge, SIGN_SQRTG
    )
    lambda_a, lambda_b = flock.calc_offset_coefficients(field.r, dr_datheta)

    nu_star_crit = flock.calc_nu_star_crit(field.r)
    lambda_s = flock.calc_finite_boundary_layer_correction(field, field.r, dr_datheta)

    print(f"s_tor: {s}")
    print(f"  Lambda_A     = {lambda_a:.6e}")
    print(f"  Lambda_B     = {lambda_b:.6e}")
    print(f"  nu_star_crit = {nu_star_crit:.6e}")
    print(f"  Lambda_S     = {lambda_s:.6e}")

    if SHOULD_CALC_SHAING_CALLEN:
        lambda_lc = flock.calc_lambda_lc(field, N_ETA, dr_datheta)
        remainder = flock.get_non_omnigenous_remainder(field, N_ETA, dr_datheta)
        print(f"  lambda_LC_bB = {lambda_lc:.6e}")
        print(f"  remainder    = {remainder:.6e}")
