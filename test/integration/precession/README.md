# Precession tests

## `test_precession_analytic`

Compares the numerical mode of the offset correction with its large-$A_m$ asymptote
for $A_1 > 200$.

$$
   A_m
   \define
   \left|{\frac{\pi R m \v_E^\ast \eta_b^2}{r \nu_\ast I_B(\eta_b)}
   \right|^{1/2}
   \left|
   2 \eta_b B^{\prime\prime}
   \right|^{-1/4}
$$

**Expected error.** The next correction to the asymptote is
$\mathcal{O}(1/\ln A_1)$, so the test uses the tolerance $0.5/\ln A_1$. The scaling
is a result of ignoring the change of the bounce time $\tau_b$ with $\eta$ near
the trapped passing boundary $\eta_b$, expressing the analytic solution only by quantities
at $\eta = \eta_b$.

**Resolution limit.** For large $l_c\hat\Omega$ the solution is confined to a thin
layer at $\eta_c$. Once the layer is narrower than the $\eta$ grid near $\eta_c$,
the result is no longer resolved and drifts away from the $1/\ln A_1$ trend.
`get_max_valid_l_c` in `precession.f90` returns the largest $l_c$ for which the
layer still contains at least two grid points. The test calls it with
$\hat\Omega = 1$ to get $\hat\Omega_\mathrm{max} = l_{c,\max}/l_c$ and only runs up
to that value.
