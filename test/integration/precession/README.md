# Precession tests

## Open issue: flux mode vs. analytic asymptote (`test_precession_analytic`)

The flux check in `test_precession_analytic` currently **fails**: the numerical
flux mode is consistently about `-0.84` times the implemented analytic estimate.

### Implemented analytic formula

Anti-sigma field `B = B_0 (1 + eps_0 cos(chi) + eps_1 cos(theta)(1 - sign cos(chi)))`
with `chi = M theta - N phi`, `eta_b = 1/(B_0 (1 + |eps_0|))`,
`B'' = |B_0 eps_0| N^2`:

```
I_B    = 4 sqrt(2|eps_0|/(1+|eps_0|)) / (N B_0^2) * (1 + 2/3 |eps_0|)
J_1    = -1/4 eps_1/eps_0 (1 + 16/3 |eps_0|) I_B B_0
A_m    = sqrt(pi m l_c Omega_hat eta_b^2 / I_B) * (2 eta_b B'')^(-1/4)
flux_1 = sqrt(2) eta_b J_1 / (A_1 sqrt(ln A_1))           (valid for A_1 >> 1)
```

### What the code solves

`compute_flux_modes` integrates (with `kappa = m l_c nfp Omega_hat / (2 (M iota - N))`)

```
d/deta (eta I_j du/deta) = i kappa tau_b u,      u regular at eta_t
flux_m = -Int_{eta_c}^{eta_t} Dr_m u/u(eta_c) deta
```

The minus sign comes from integrating from deeply (`t = sqrt(eta_t - eta_c)`) to
barely trapped (`t = 0`). It is the convention that `test_g_zero_precession`
checks against the fieldline integrals.

What has been verified:

- An independent Python solve of the same equation from the dumped bounce data
  reproduces the Fortran flux to `<= 3e-4` relative, for `Omega_hat` = 2 ... 1.4e9.
- The inputs match the analytic ingredients: `I_j(eta_c) = I_B` (0.6339 vs
  0.6338), and `tau_b` follows the barely trapped log up to a constant shift.
  The drift is `Dr_1(eta_c) = 4/3 J_1`, where `4/3 = (B_phi + iota B_theta)/(0.75 psi_edge)`
  is the conversion factor of `test_g_zero_precession`.

### Asymptote of the code's equation

Near `eta_c` the coefficients are `eta I_j ~ eta_b I_B` and `tau_b ~ -a ln(x/eta_b)`,
with `x = eta - eta_c` and `a = 4 eta_b (2 eta_b B'')^(-1/2)`. This gives a
boundary layer `u/u(eta_c) = exp(-sqrt(i Q L) x)`, where `L = ln(1/x) ~ ln A`
and `Q = 2 A^2/(pi eta_b^2)`. It follows that `Re Int u/u(eta_c) dx = 1/sqrt(2 Q L)`, so

```
flux_1 -> -Dr_1(eta_c) eta_b sqrt(pi) / (2 A_1 sqrt(ln A_1))
flux_1 / flux_1(implemented) -> -sqrt(pi/8) * Dr_1/J_1 = -sqrt(pi/8) * 4/3 = -0.836
```

An idealized solve (constant `I_B`, pure log `tau_b`, constant drift) confirms
this limit. The approach is slow, `O(1/ln A)`. The measured ratios
(`eps_0 = -0.0125`, `eps_1 = -0.0005`, `N = 1`) are:

| A_1     | numeric / implemented |
|---------|-----------------------|
| 2.1e2   | -0.687                |
| 1.0e3   | -0.784                |
| 5.6e3   | -0.800                |
| 3.0e4   | -0.804                |
| 5.8e5   | -0.816                |

### Suspected origin of the factor `-sqrt(pi/8) * 4/3`

1. **Sign:** the code uses `flux = -Int_{eta_c}^{eta_t} ...`, which is fixed by
   `test_g_zero_precession`. The formula appears to assume `+Int`, or a `J_1`
   with the opposite sign.
2. **Factor `3/4`:** `J_1` seems to be in `I_V` units, without the conversion factor
   `(B_phi + iota B_theta)/(0.75 psi_edge)` that the code's radial drift carries.
   The code has `Dr_1 = 4/3 J_1`.
3. **Factor `sqrt(8/pi)`:** this is most likely in the derivation of the
   boundary-layer integral or in the `pi` inside `A_m`. The exponential layer gives
   `sqrt(pi)/2`, where the formula has `sqrt(2)`. Because the factor contains `pi`,
   it cannot come from a simple prefactor in the code (`1/2` in `kappa`, `nfp`,
   half-bounce `I_j`, `0.75`).

Fixing (1) and (2) alone would leave a ratio of `sqrt(pi/8) = 0.627`.
