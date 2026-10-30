# Precession tests

## Open issue: flux mode vs. analytic asymptote (`test_precession_analytic`)

The flux check in `test_precession_analytic` currently **fails**: the numerical
flux mode is consistently about $-0.84$ times the implemented analytic estimate.

### Implemented analytic formula

The test uses the anti-sigma field

```math
B = B_0\left(1 + \varepsilon_0\cos\chi + \varepsilon_1\cos\vartheta\,(1 - \sigma\cos\chi)\right),
\qquad \chi = M\vartheta - N\varphi,
```

with $\eta_b = 1/\big(B_0(1+|\varepsilon_0|)\big)$ and $B'' = |B_0\varepsilon_0| N^2$:

```math
\begin{aligned}
I_B &= \frac{4}{N B_0^2}\sqrt{\frac{2|\varepsilon_0|}{1+|\varepsilon_0|}}
       \left(1 + \tfrac{2}{3}|\varepsilon_0|\right), \\
J_1 &= -\frac{1}{4}\frac{\varepsilon_1}{\varepsilon_0}
       \left(1 + \tfrac{16}{3}|\varepsilon_0|\right) I_B B_0, \\
A_m &= \sqrt{\frac{\pi\, m\, l_c\, \hat\Omega\, \eta_b^2}{I_B}}\;
       \left(2\eta_b B''\right)^{-1/4}, \\
\Gamma_1^\mathrm{impl} &= \frac{\sqrt{2}\,\eta_b J_1}{A_1\sqrt{\ln A_1}}
       \qquad (A_1 \gg 1).
\end{aligned}
```

### What the code solves

With $\kappa = \dfrac{m\, l_c\, n_\mathrm{fp}\, \hat\Omega}{2(M\iota - N)}$,
`compute_flux_modes` integrates

```math
\frac{d}{d\eta}\left(\eta I_j \frac{du}{d\eta}\right) = i\kappa\,\tau_b\, u,
\qquad u \text{ regular at } \eta_t,
```

```math
\Gamma_m = -\int_{\eta_c}^{\eta_t} \Delta r_m\, \frac{u}{u(\eta_c)}\, d\eta .
```

The minus sign comes from integrating from deeply trapped
($t = \sqrt{\eta_t - \eta_c}$) to barely trapped ($t = 0$). It is the convention
that `test_g_zero_precession` checks against the fieldline integrals.

What has been verified:

- An independent Python solve of the same equation from the dumped bounce data
  reproduces the Fortran flux to $\le 3\cdot 10^{-4}$ relative, for
  $\hat\Omega = 2 \dots 1.4\cdot 10^{9}$.
- The inputs match the analytic ingredients: $I_j(\eta_c) = I_B$ (0.6339 vs
  0.6338), and $\tau_b$ follows the barely trapped log up to a constant shift.
  The drift is $\Delta r_1(\eta_c) = \tfrac{4}{3} J_1$, where
  $\tfrac{4}{3} = (B_\varphi + \iota B_\vartheta)/(0.75\,\psi_\mathrm{edge})$
  is the conversion factor of `test_g_zero_precession`.

### Asymptote of the code's equation

Near $\eta_c$, with $x = \eta - \eta_c$, the coefficients are
$\eta I_j \approx \eta_b I_B$ and $\tau_b \approx -a\ln(x/\eta_b)$, where
$a = 4\eta_b(2\eta_b B'')^{-1/2}$. This gives a boundary layer

```math
\frac{u}{u(\eta_c)} \approx e^{-\sqrt{iQL}\,x},
\qquad L = \ln\frac{1}{x} \approx \ln A,
\qquad Q = \frac{2A^2}{\pi\eta_b^2},
```

so that $\mathrm{Re}\int u/u(\eta_c)\,dx = 1/\sqrt{2QL}$ and

```math
\Gamma_1 \to -\Delta r_1(\eta_c)\,\frac{\eta_b\sqrt{\pi}}{2A_1\sqrt{\ln A_1}},
\qquad
\frac{\Gamma_1}{\Gamma_1^\mathrm{impl}} \to
-\sqrt{\frac{\pi}{8}}\,\frac{\Delta r_1}{J_1}
= -\sqrt{\frac{\pi}{8}}\cdot\frac{4}{3} \approx -0.836 .
```

An idealized solve (constant $I_B$, pure log $\tau_b$, constant drift) confirms
this limit. The approach is slow, $\mathcal{O}(1/\ln A)$. The measured ratios
($\varepsilon_0 = -0.0125$, $\varepsilon_1 = -0.0005$, $N = 1$) are:

| $A_1$            | $\Gamma_1/\Gamma_1^\mathrm{impl}$ |
|------------------|-----------------------------------|
| $2.1\cdot10^{2}$ | $-0.687$                          |
| $1.0\cdot10^{3}$ | $-0.784$                          |
| $5.6\cdot10^{3}$ | $-0.800$                          |
| $3.0\cdot10^{4}$ | $-0.804$                          |
| $5.8\cdot10^{5}$ | $-0.816$                          |

### Suspected origin of the factor $-\sqrt{\pi/8}\cdot 4/3$

1. **Sign:** the code uses $\Gamma = -\int_{\eta_c}^{\eta_t}\dots$, which is fixed
   by `test_g_zero_precession`. The formula appears to assume $+\int$, or a
   $J_1$ with the opposite sign.
2. **Factor $3/4$:** $J_1$ seems to be in $I_\mathcal{V}$ units, without the
   conversion factor $(B_\varphi + \iota B_\vartheta)/(0.75\,\psi_\mathrm{edge})$
   that the code's radial drift carries. The code has $\Delta r_1 = \tfrac{4}{3}J_1$.
3. **Factor $\sqrt{8/\pi}$:** this is most likely in the derivation of the
   boundary-layer integral or in the $\pi$ inside $A_m$. The exponential layer
   gives $\sqrt{\pi}/2$, where the formula has $\sqrt{2}$. Because the factor
   contains $\pi$, it cannot come from a simple prefactor in the code
   ($\tfrac12$ in $\kappa$, $n_\mathrm{fp}$, half-bounce $I_j$, $0.75$).

Fixing (1) and (2) alone would leave a ratio of $\sqrt{\pi/8} \approx 0.627$.
