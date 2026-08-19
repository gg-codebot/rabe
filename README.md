<img width="0" height="20" align="right">
<img src="assets/rabe_logo.png#gh-light-mode-only" alt="rabe logo by Georg Grassler" width="200" align="right"><img src="assets/rabe_logo_dark.png#gh-dark-mode-only" alt="rabe logo by Georg Grassler" width="200" align="right">

# rabe

[![Documentation](https://img.shields.io/badge/docs-FORD-blue)](https://itpplasma.github.io/rabe/)

An implementation of the nea**r**-omnigenous, **a**symptotic **b**ootstrap **e**xpressions of Ref. [1,2].

## Summary

In case of the Lorentz collision model,
the mono-energetic bootstrap and Ware pinch coefficients, $\bar D_{31}$
and $\bar D_{13}$, can be expressed
via a dimensionless geometry parameter $\lambda_{bB}$,

$$
\bar D_{31}=\bar D_{13} = \frac{1}{3} v \rho_L B \lambda_{bB}, \tag{1}
$$

scaling out particle velocity $v$, $\rho_L$ and $B$. This parameter can
be split into a ''symmetric'' contribution and a ''non-symmetric''  off-set,

$$
\lambda_{bB} = \lambda_\mathrm{sym} + \lambda_\mathrm{off}. \tag{2}
$$

For near-omnigenous stellarators at low plasma collisionality
$\nu_\ast=\pi R \nu_{\perp}/v$, where $\nu_\perp$ and $R$ are deflection frequency and
device major radius, respectively, the off-set can be expressed as

$$
    \lambda_{\text{off}} = \frac{\Lambda_\mathrm{A}}{\sqrt{\nu_\ast}} + \frac{\Lambda_\mathrm{B}}{\nu_\ast}, \tag{3}
$$

where $\Lambda_\mathrm{A}$ and $\Lambda_\mathrm{B}$ are the geometrical factors
due to the variation of the trapped-passing boundary layer width and the misalignment
of local maxima, respectively (see Ref. [1,2]).

Given a VMEC equilibrium file, `rabe` outputs those geometric coefficients
needed to evaluate $\lambda_\mathrm{off}$ at any collisionality.
Additionally, it can compute also the coefficients for the
''symmetric'' part

$$
\lambda_\mathrm{sym} = \lambda_{bB}^\mathrm{SC} + \lambda_{bB}^\mathrm{HGM}, \qquad \lambda_{bB}^\mathrm{HGM} = \Lambda_\mathrm{S}\sqrt{\nu_\ast} \tag{4}
$$

where the Shaing-Callen asymptotic $\lambda_{bB}^\mathrm{SC}$ is
approximated by the perfect omnigenous asymptotic
$\lambda_{bB}^\mathrm{SC} \rightarrow \lambda_{bB}^\mathrm{LC}$
obtained by Landreman and Catto [3] and
$\lambda_{bB}^\mathrm{HGM}$ is the contribution due to the finite
boundary correction derived by Helander, Geiger and Maasberg [4].

## Prerequisites

- Fortran compiler (gfortran)
- CMake >= 3.24

On Debian/Ubuntu:

```bash
sudo apt-get install gfortran cmake
```

On macOS (Homebrew):

```bash
brew install gfortran
export FC=$(ls $(brew --prefix)/bin/gfortran-* | sort -V | tail -1)
```

Make sure that `FC` points to your installed gfortran compiler.

## Example

This example walks through a complete run using the QH stellarator equilibrium
from Ref. [5], which ships with the repository as a test input.
Commands are given for `bash` shell. While in `rabe`
run

**Step 1 — build:**

```bash
make clean
make CONFIG=Release
```

This builds the executable `rabe.x` in Release mode and writes it to `build`.

**Step 2 — create a working directory and link the inputs:**

```bash
mkdir run_example && cd run_example
ln -s ../test/golden/input/rabe.in .
ln -s ../test/integration/vmec/input/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc .
```

For explanation of the input parameters see the Input / Output section below.

**Step 3 — run:**

```bash
../build/rabe.x
```

`rabe` reads `rabe.in` and `field_file` from the working directory and writes
results to `rabe.nc` (NetCDF) and `rabe.dat` (plain text) after completion.

**Step 4 — visualize:**

```bash
python ../test/golden/plot_golden.py rabe.nc
```

or with Octave:

```bash
octave ../test/golden/plot_golden.m rabe.nc
```

Both produce `rabe_output.png` showing the off-set factors ($\Lambda_A$,
$\Lambda_B$), $\Lambda_\mathrm{S}$, and — if `should_calc_shaing_callen =
.true.` — the Shaing-Callen asymptotic, all plotted against $s_\mathrm{tor}$.
The Python script requires `matplotlib` and `xarray`. The Octave script
requires the `netcdf` package (`pkg install -forge netcdf` if not present).

## Python interface

In addition to the Fortran executable, `rabe` ships a Python package with the
same name that exposes the core computation through a thin f90wrap binding.

### Installation

As a wheel (no compiler required if a matching binary is available):

```bash
pip install rabe
```

Or directly from source (requires the same Fortran prerequisites
listed under Prerequisites above):

```bash
pip install -e .
```

### Examples

Ready-to-run scripts are in the `python/` directory:

| Script | Field type | Use case |
| --- | --- | --- |
| `python/example.py` | `BoozerField` | Real VMEC equilibrium from a `.nc` file |
| `python/example_fourier.py` | `FourierField` | Analytical Fourier-mode field |

**`BoozerField`** (`boozer_field_t`): loads a VMEC NetCDF 3D equilibrium and
converts it to Boozer coordinates. This is the same field representation used
by the executable.

**`FourierField`** (`fourier_field_t`): builds an 2D field directly from
a flat list of Boozer Fourier modes plus surface values for covariant components.
Useful for benchmark and optimisation studies.

Both field types are accepted by the main downstream API via `FlockOfFieldlines`:

```python
from rabe.fieldline_mod import FlockOfFieldlines
flock = FlockOfFieldlines(max_n_fieldlines, iota, field, M_pol, N_tor, nfp)
lambda_a, lambda_b = flock.calc_offset_coefficients(R, dr_dAtheta)
nu_star_crit       = flock.calc_nu_star_crit(R)
```

## Input / Output

`rabe` reads its configuration from a namelist file `rabe.in` in the working directory:

```fortran
&rabe_config
    field_file = "wout_example.nc",      ! VMEC equilibrium file (NetCDF)
    M_pol = -1.0,                        ! poloidal helicity of omnigenity
    N_tor = 4.0,                         ! toroidal helicity of omnigenity
    s_tor = 0.25, 0.5, 0.75,             ! explicit surface list OR
!   s_tor_min = 0.1,                     ! \
!   s_tor_max = 0.9,                     !  } n_s_tor equi-spaced surfaces
!   n_s_tor   = 9,                       ! /
    sign_sqrtg = -1.0,                   ! sign of the Jacobian sqrt(g)
    max_n_fieldlines = 200,              ! maximum field lines per surface
    should_calc_shaing_callen = .true.,  ! compute $\lambda_{bB}^\mathrm{LC}$
    n_eta = 100,                         ! level resolution for trapped
                                         ! particle fraction computation
    unsafe_mode = .false.                ! if .true., NaN-fill outputs for
                                         ! surfaces that fail a sanity check
                                         ! instead of aborting execution
/
```

The type of omnigenity is set by the helicity

| type of omnigenity | `M_pol` | `N_tor` |
| --- | --- | --- |
| quasi-axisymmetric | 1 | 0 |
| quasi-isodynamic | 0 | $N_p$ |
| quasi-helicalsymmetric | $M$ | $N_p$ |

where $N_p$ is the number of field periods and $M$ is
given by the poloidal mode number $m$ of the strongest mode
$(m,n=N_p)$ of the Boozer spectrum

$$
 B(\vartheta, \varphi) = \sum_{m,n} B_{m,n} \cos{(m\vartheta - n\varphi)}. \tag{4}
$$

The list of surfaces on which to compute can either be given explicitly in `s_tor` OR via a uniform range. `sign_sqrtg` is globally applied
to all output coefficients to account for different coordinate conventions.
It should be set to the same value as `signgs` in the `VMEC` output `wout_*.nc` and is usually `sign_sqrtg=-1.0`.

By default (`unsafe_mode = .false.`), any failed sanity check halts the run
immediately with an error. This is the recommended behaviour as those checks point out not suited inputs e.g. violation of stellarator symmetry.
Setting `unsafe_mode = .true.` allows the run to continue, but outputs
are set to `NaN` for any surface where a check failed.

Results are written to `rabe.nc` (NetCDF) and `rabe.dat` (plain text), with one
value per flux surface. Both files contain the same variables:

| Variable | Description |
| --- | --- |
| `s_tor` | toroidal flux $\psi$ normalized to flux at edge $\psi_a$ as $s_\mathrm{tor} = \psi/\psi_a$ |
| `Lambda_A` | $1/\sqrt{\nu_\ast}$ factor ($\Lambda_\mathrm{A}$ in Eq. 3) |
| `Lambda_B` | $1/\nu_\ast$ factor ($\Lambda_\mathrm{B}$ in Eq. 3) |
| `nu_star_crit` | lower collisionality limit for asymptotic model validity |
| `Lambda_S` | $\sqrt{\nu_\ast}$ correction for finite boundary layer width |
| `split_maxima` | 1 if omnigeneity violation is too strong, 0 otherwise |
| `R` | major radius [m] (reference length scale for $\nu_\ast$) |

`split_maxima` warns the user to treat the results with caution. If the
violation of omnigenity is too strong, local maxima contours are not merely
deformed, but also get split, which this flag notes.

For fast runs suited for **optimization**, we recommend to set
`should_calc_shaing_callen = .false.`. If it is enabled, the output
is extended by

| Variable | Description |
| --- | --- |
| `lambda_LC_bB` | omnigenous Shaing-Callen coefficient $\lambda_{bB}^\mathrm{LC}$ |
| `remainder` | non-omnigenous remainder (only as a prototype) |

The `remainder` is a proxy of how much $\lambda_{bB}^\mathrm{LC}$ differs
from the actual $\lambda_{bB}^\mathrm{SC}$, but it does not include the
effect of bootstrap resonances and is not yet fully validated.

## Adjoint Gradients

**Status: under active development.**

Adjoint-based gradients of the off-set coefficients with respect
to e.g. Boozer modes $B_{mn}$ and $\iota$ are being implemented internally currently.
A usage guide will be added here once the feature lands.

## Build and Tests

We use CMake for build configuration. If it is available on your machine, we
recommend [Ninja](https://ninja-build.org) as the generator

```bash
export CMAKE_GENERATOR=Ninja
```

Run

```bash
make
```

to build the executable (debug build) or `make CONFIG=Release` for an optimized
build. You can run the main suite of tests used to ensure correctness
of the code with

```bash
make test
```

or run all tests (including ones that take quite some time)

```bash
make test_all
```

The tests as well as their description, can be found in `test`.

### Overriding the libneo dependency

By default the build fetches a pinned libneo commit. Two explicit options change this:

- `-DLIBNEO_REF=<branch|tag|sha>` selects a different git ref:

  ```bash
  make LIBNEO_REF=main
  # or directly: cmake -S . -B build -DLIBNEO_REF=main
  ```

- `-DLIBNEO_PATH=<dir>` uses a local checkout instead of fetching:

  ```bash
  make LIBNEO_PATH=/path/to/libneo
  # or directly: cmake -S . -B build -DLIBNEO_PATH=/path/to/libneo
  ```

## Third Party

Fetched automatically during build:

- [libneo](https://github.com/itpplasma/libneo) for field file I/O; pass `-DLIBNEO_PATH=<dir>` to use a local checkout instead (MIT)
- [`quadpack`](https://github.com/jacobwilliams/quadpack) for numerical integration (BSD-3-Clause)
- [`fortio`](https://github.com/lazy-fortran/fortio) for NetCDF input/output (MIT)
- [`pyplot-fortran`](https://github.com/jacobwilliams/pyplot-fortran) optional for visualization; source and license under `plot_lib` (BSD-3-Clause)

## Citing

If you use this code in any of your studies, please cite Refs [1,2].

See also [`CITATION.cff`](CITATION.cff) for machine-readable citation metadata.

## References

[1] Georg S. Grassler, Christopher G. Albert, Sergei V. Kasilov, Winfried Kernbichler, *Asymptotic modeling of the bootstrap and Ware pinch effect in near omnigenous stellarators*, Phys. Plasmas; 33 (7): 070703. [doi:10.1063/5.0332431](https://doi.org/10.1063/5.0332431) (2026)

[2] C.G Albert et al., *On the convergence of bootstrap current to the Shaing–Callen limit in stellarators*, Journal of Plasma Physics, 91(3), p. E77. [doi:10.1017/S0022377825000200](https://doi.org/10.1017/S0022377825000200) (2025)

[3] M. Landreman and P. J. Catto, *Omnigenity as generalized quasisymmetry*, Phys. Plasmas 19, 056103 [doi.org/10.1063/1.3693187](https://doi.org/10.1063/1.3693187) (2012)

[4] P. Helander, J. Geiger, and H. Maassberg, “On the bootstrap current in
stellarators and tokamaks”, Phys. Plasmas 18, 092505 [doi.org/10.1063/1.3633940](https://doi.org/10.1063/1.3633940) (2011)

[5] M. Landreman et al., *Optimization of quasi-symmetric stellarators with self-consistent bootstrap current and energetic particle confinement*, Phys. Plasmas 29, [doi:10.1063/5.0098166](https://doi.org/10.1063/5.0098166) (2022)

[6] John R. Cary & Svetlana G. Shasharina, *Omnigenity and quasihelicity in helical plasma confinement systems, Phys. Plasmas 4*, 3323–3333, [doi:10.1063/1.872473](https://doi.org/10.1063/1.872473) (1997)
