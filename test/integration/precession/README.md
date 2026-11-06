# Precession tests

## `test_precession_analytic`

Compares the numerical flux mode $\Gamma_1$ with its large-$A$ asymptote
$\Gamma_1^\mathrm{an}$ for $A_1 > 200$.

**Expected error.** The next correction to the asymptote is
$\mathcal{O}(1/\ln A_1)$, so the test uses the tolerance $0.5/\ln A_1$. The
measured ratio $\Gamma_1/\Gamma_1^\mathrm{an}$ is about $0.95$ and approaches $1$
as $A_1$ grows.

**Resolution limit.** For large $l_c\hat\Omega$ the solution is confined to a thin
layer at $\eta_c$. Once the layer is narrower than the $\eta$ grid near $\eta_c$,
the result is no longer resolved and drifts away from the $1/\ln A_1$ trend.
`get_max_valid_l_c` in `precession.f90` returns the largest $l_c$ for which the
layer still contains at least two grid points. The test calls it with
$\hat\Omega = 1$ to get $\hat\Omega_\mathrm{max} = l_{c,\max}/l_c$ and only runs up
to that value.
