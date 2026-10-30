# Precession tests

## Open issue: flux mode vs. analytic asymptote (`test_precession_analytic`)

The flux check in `test_precession_analytic` currently **fails**: the numerical
flux mode is consistently about $-1.22$ times the implemented analytic estimate.
The asymptotic limit is $-\sqrt{\pi/2} \approx -1.25$.

### Implemented analytic formula

The test uses the anti-sigma field

```math
B = B_0\left(1 + \varepsilon_0\cos\chi + \varepsilon_1\cos\vartheta\,(1 - \sigma\cos\chi)\right),
\qquad \chi = M\vartheta - N\varphi,
```

with $\iota = 0$, so that $\vartheta$ is constant along a well. It uses
$\eta_b = 1/\big(B_0(1+|\varepsilon_0|)\big)$ and $B'' = |B_0\varepsilon_0| N^2$:

```math
\begin{aligned}
I_B &= \frac{4}{N B_0^2}\sqrt{\frac{2|\varepsilon_0|}{1+|\varepsilon_0|}}
       \left(1 + \tfrac{2}{3}|\varepsilon_0|\right), \\
J_1 &= -\frac{2}{N B_0}\sqrt{\frac{2|\varepsilon_0|}{1+|\varepsilon_0|}}\,
       \frac{\varepsilon_1}{\varepsilon_0}\left(1 + 3|\varepsilon_0|\right), \\
A_m &= \sqrt{\frac{\pi\, m\, l_c\, \hat\Omega\, \eta_b^2}{I_B}}\;
       \left(2\eta_b B''\right)^{-1/4}, \\
\Gamma_1^\mathrm{impl} &= \frac{\sqrt{2}\,\eta_b J_1}{A_1\sqrt{\ln A_1}}
       \qquad (A_1 \gg 1).
\end{aligned}
```

Here $I_B = \langle I_j(\eta_c)\rangle_{\xi_0}$ is the field line average of $I_j$
at the barely trapped boundary. It is independent of $\iota$ to
$\mathcal{O}(\varepsilon_1^2)$. $J_1 = \tfrac{2}{3}\,\partial_\eta I_{\mathcal{V},1}$
is the $m=1$ mode at $\eta_c$. The barely trapped
$\partial_\eta I_{\mathcal{V},1}$ for $\iota \approx 0$ is verified in
`test_radial_drift`. Unlike $I_B$, $J_1$ depends on $\iota$: at $\iota = 1$ the
drift mode is about $2.9\times$ smaller.

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
- The inputs match the analytic ingredients:
  - $I_j(\eta_c) = I_B$ (0.63381 vs 0.63378).
  - $\tau_b$ follows the barely trapped log up to a constant shift.
  - $\Delta r_1(\eta_c) = 2 J_1$: $-0.02633$ vs $2\cdot(-0.01304)$, consistent with
    $\Delta r = \partial_\eta I_\mathcal{V}/0.75$ for $B_\varphi + \iota B_\vartheta = 1$
    and $\psi_\mathrm{edge} = 1$.

### Asymptote of the code's equation

Near $\eta_c$, with $x = \eta - \eta_c$, the coefficients are
$\eta I_j \approx \eta_b I_B$ and $\tau_b \approx -a\ln(x/\eta_b)$, where
$a = 4\eta_b(2\eta_b B'')^{-1/2}$. This gives a boundary layer

```math
\frac{u}{u(\eta_c)} \approx e^{-\sqrt{iQL}\,x},
\qquad L = \ln\frac{1}{x} \approx \ln A,
\qquad Q = \frac{2A^2}{\pi\eta_b^2},
```

so that $\mathrm{Re}\int u/u(\eta_c)\,dx = 1/\sqrt{2QL}$ and, with $\Delta r_1 = 2J_1$,

```math
\Gamma_1 \to -\Delta r_1(\eta_c)\,\frac{\eta_b\sqrt{\pi}}{2A_1\sqrt{\ln A_1}}
= -\frac{\sqrt{\pi}\,\eta_b J_1}{A_1\sqrt{\ln A_1}},
\qquad
\frac{\Gamma_1}{\Gamma_1^\mathrm{impl}} \to -\sqrt{\frac{\pi}{2}} \approx -1.253 .
```

An idealized solve (constant $I_B$, pure log $\tau_b$, constant drift) confirms
this limit. The approach is slow, $\mathcal{O}(1/\ln A)$. The measured ratios
($\varepsilon_0 = -0.0125$, $\varepsilon_1 = -0.0005$, $N = 1$, $\iota = 0$) are:

| $A_1$            | $\Gamma_1/\Gamma_1^\mathrm{impl}$ |
|------------------|-----------------------------------|
| $2.1\cdot10^{2}$ | $-1.223$                          |
| $1.8\cdot10^{3}$ | $-1.208$                          |
| $1.7\cdot10^{4}$ | $-1.215$                          |
| $1.7\cdot10^{5}$ | $-1.219$                          |
| $5.8\cdot10^{5}$ | $-1.233$                          |

### Remaining discrepancy: $-\sqrt{\pi/2}$

1. **Sign:** the code uses $\Gamma = -\int_{\eta_c}^{\eta_t}\dots$, which is fixed
   by `test_g_zero_precession`. The formula appears to assume $+\int$.
2. **Factor $\sqrt{\pi/2}$:** the boundary-layer integral of the code's equation
   gives $\sqrt{\pi}\,\eta_b J_1$, where the implemented formula has
   $\sqrt{2}\,\eta_b J_1$. The suspects are the derivation of that integral and
   the $\pi$ inside $A_m$.

### Resolved along the way

- The previous $J_1 = -\tfrac14\frac{\varepsilon_1}{\varepsilon_0}(1+\tfrac{16}{3}|\varepsilon_0|)\,I_B B_0$
  was half of the verified $\iota = 0$ estimate. The test also ran at $\iota = 1$,
  where $J_1$ does not apply.
- `compute_bounce_integrals` estimated the deeply trapped bounce time with
  $h^\vartheta$ ($= \iota h^\varphi$) instead of $h^\varphi$. This gave an
  infinite time step at $\iota = 0$, and a step $\sim 100\times$ too large at
  $\iota = 0.01$.
- The flux ODE overflowed to NaN for $A \gtrsim 1.5\cdot 10^4$. It now
  integrates in renormalized segments.
