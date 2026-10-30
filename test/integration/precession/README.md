# Precession tests

## Flux mode vs. analytic asymptote (`test_precession_analytic`)

`test_precession_analytic` checks the numerical flux mode against the large-$A$
asymptote for $A_1 > 200$, with tolerance $1/\ln A_1$ (the order of the next
correction). The measured ratio is $0.95 \dots 0.96$ and approaches $1$ as
$\mathcal{O}(1/\ln A)$.

### Analytic formula

The test uses the anti-sigma field

```math
B = B_0\left(1 + \varepsilon_0\cos\chi + \varepsilon_1\cos\vartheta\,(1 - \sigma\cos\chi)\right),
\qquad \chi = M\vartheta - N\varphi,
```

with $\iota = 0$, so that $\vartheta$ is constant along a well. It uses
$\eta_b = 1/\big(B_0(1+|\varepsilon_0|)\big)$, $B'' = |B_0\varepsilon_0| N^2$, and
$l_c = \pi R/(2\nu^*)$, where the $\pi$ is part of $l_c$:

```math
\begin{aligned}
I_B &= \frac{4}{N B_0^2}\sqrt{\frac{2|\varepsilon_0|}{1+|\varepsilon_0|}}
       \left(1 + \tfrac{2}{3}|\varepsilon_0|\right), \\
J_1 &= \frac{2}{N B_0}\sqrt{\frac{2|\varepsilon_0|}{1+|\varepsilon_0|}}\,
       \frac{\varepsilon_1}{\varepsilon_0}\left(1 + 3|\varepsilon_0|\right), \\
A_m &= \sqrt{\frac{2\, m\, l_c\, \hat\Omega\, \eta_b^2}{I_B}}\;
       \left(2\eta_b B''\right)^{-1/4}, \\
\Gamma_1^\mathrm{an} &= \frac{\sqrt{2}\,\eta_b J_1}{A_1\sqrt{\ln A_1}}
       \qquad (A_1 \gg 1).
\end{aligned}
```

Here $I_B = \langle I_j(\eta_c)\rangle_{\xi_0}$ is the field line average of $I_j$
at the barely trapped boundary. It is independent of $\iota$ to
$\mathcal{O}(\varepsilon_1^2)$ and is checked against the code.
$J_1 = -\tfrac{2}{3}\,\partial_\eta I_{\mathcal{V},1}$ is the $m=1$ mode at $\eta_c$.
The barely trapped $\partial_\eta I_{\mathcal{V},1}$ for $\iota \approx 0$ is
verified in `test_radial_drift`. Unlike $I_B$, $J_1$ depends on $\iota$: at
$\iota = 1$ the drift mode is about $2.9\times$ smaller.

### What the code solves

With $\kappa = \dfrac{m\, l_c\, n_\mathrm{fp}\, \hat\Omega}{2(M\iota - N)}$,
`compute_flux_modes` integrates

```math
\frac{d}{d\eta}\left(\eta I_j \frac{du}{d\eta}\right) = i\kappa\,\tau_b\, u,
\qquad u \text{ regular at } \eta_t,
```

```math
\Gamma_m = -\int_{\eta_c}^{\eta_t} \Delta s_m\, \frac{u}{u(\eta_c)}\, d\eta ,
```

where $\Delta s_m$ is the $m$-th mode (over field line labels) of the change in
$s$ over one bounce, normalized by the Larmor radius (`radial_drift` in the
code).

The minus sign comes from integrating from deeply trapped
($t = \sqrt{\eta_t - \eta_c}$) to barely trapped ($t = 0$). It is the convention
that `test_g_zero_precession` checks against the fieldline integrals.

What has been verified:

- An independent Python solve of the same equation from the dumped bounce data
  reproduces the Fortran flux to $\le 3\cdot 10^{-4}$ relative.
- The inputs match the analytic ingredients:
  - $I_j(\eta_c) = I_B$ (0.63381 vs 0.63378).
  - $\tau_b$ follows the barely trapped log up to a constant shift.
  - $\Delta s_1(\eta_c) = -2 J_1$: $-0.02633$ vs $-2\cdot 0.01304$.

### Asymptote of the code's equation

Near $\eta_c$, with $x = \eta - \eta_c$, the coefficients are
$\eta I_j \approx \eta_b I_B$ and $\tau_b \approx -a\ln(x/\eta_b)$, where
$a = 4\eta_b(2\eta_b B'')^{-1/2}$. This gives a boundary layer
$u/u(\eta_c) \approx e^{-\sqrt{iQL}\,x}$ with $L \approx \ln A$ and
$Q = |\kappa| a/(\eta_b I_B) = A^2/\eta_b^2$, so that

```math
\Gamma_1 \to -\Delta s_1(\eta_c)\,\frac{\eta_b}{\sqrt{2}\,A_1\sqrt{\ln A_1}}
= \frac{\sqrt{2}\,\eta_b J_1}{A_1\sqrt{\ln A_1}} = \Gamma_1^\mathrm{an}.
```

The measured ratios ($\varepsilon_0 = -0.0125$, $\varepsilon_1 = -0.0005$,
$N = 1$, $\iota = 0$) are:

| $A_1$            | $\Gamma_1/\Gamma_1^\mathrm{an}$ |
|------------------|---------------------------------|
| $2.1\cdot10^{2}$ | $0.951$                         |
| $1.0\cdot10^{3}$ | $0.948$                         |
| $5.6\cdot10^{3}$ | $0.955$                         |
| $3.0\cdot10^{4}$ | $0.960$                         |
| $1.7\cdot10^{5}$ | $0.964$                         |

### Resolved along the way

- $J_1$ had the wrong sign: $J_1 = -\tfrac23\,\partial_\eta I_{\mathcal{V},1} > 0$.
- $A_m$ lacked the factor 2 and had a spurious $\pi$, which is already
  contained in $l_c$.
- The previous $J_1 = -\tfrac14\frac{\varepsilon_1}{\varepsilon_0}(1+\tfrac{16}{3}|\varepsilon_0|)\,I_B B_0$
  was half of the verified $\iota = 0$ estimate. The test also ran at $\iota = 1$,
  where $J_1$ does not apply.
- `compute_bounce_integrals` estimated the deeply trapped bounce time with
  $h^\vartheta$ ($= \iota h^\varphi$) instead of $h^\varphi$. This gave an
  infinite time step at $\iota = 0$, and a step $\sim 100\times$ too large at
  $\iota = 0.01$.
- The flux ODE overflowed to NaN for $A \gtrsim 1.5\cdot 10^4$. It now
  integrates in renormalized segments.
