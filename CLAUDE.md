# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## General

- do not change the rules in this section, in `Contributing`, `Documentation`
- rules in the above mentioned sections take precedent
- if you wanna add rules yourself, that can be changed, make a new section
- those new rules are always superseeded if applicable by the rules in the above sections
- if you think a rule is worth to be updated or to be added to those sections, suggest it in the description of an pr.
- this repo is inteded to be read and used by humans: less is more - if one writes e.g. a test, it should be impactful and docs should be short and to the point

# Contributing

- Work on feature branches, test your work regularly, commit in small chunks.
- Open a merge request for review and tag the maintainer — do not commit directly to `main`.
- maintainer tag: @GeorgGrassler

# Documentation

- be minimalistic in documentation, let layout and code speak for themself
- do not write information that can be outdated fast like specific file or variable names that are physically far removed from their mention e.g. in app about module implementation details

## Project Overview

RABE (Redl's Analytical Bootstrap Estimation) is a scientific computing project for analytical bootstrap current estimation in stellarators, implemented in modern Fortran with Python utilities for analysis and visualization.

## Build and Test Commands

### Building
- `make` - Build the project (Debug mode by default)
- `make CONFIG=Release` - Build in Release mode
- `make clean` - Clean build directory
- `make install` - Build and install

### Testing
- `make test` - Run quick tests (label: `quick`)
- `make test_slow` - Run performance-intensive tests (label: `slow`)
- `make test_all` - Run quick + slow tests
- `make plot` - Run plotting/visualization tests (label: `plot`)
- `make current` - Run tests marked as current development focus (label: `current`)
- `make internal` - Run private/internal tests (label: `internal`; requires `test/internal` submodule)
- `make golden` - Run golden record tests (label: `golden`)
- `make golden_run` - Run golden record execution only (no comparison; used for timing)
- `make golden_update` - Update golden record expected output

Single test: `ctest -R test_name --test-dir build`

## Architecture Overview

### Libraries (`src/`)

```
rabe_lib (main static library)
├── utils_lib (src/utils/)
│   constants.f90, utils.f90, diophantine.f90, find_extrema.f90,
│   fourier.f90, integrate.f90
├── field_lib (src/field/)
│   field_base.f90    - Abstract type field_t (deferred: compute_B_mod,
│                       compute_nabla_s, compute_B_sqrtg_dB_dx, compute_B_and_dB_dx,
│                       get_covariant_components)
├── vmec_lib (src/vmec/)
│   boozer_field.f90  - Concrete field_t type boozer_field_t with one loader
│                       per input format: init_from_vmec (VMEC wout),
│                       init_from_boozmn (boozmn), init_from_chartmap
│                       (extended chartmap); loaders live in libneo
│                       (boozer_sub, boozer_chartmap, boozmn_reader)
├── fieldline_lib (src/fieldline/) — see src/fieldline/README.md
│   fieldline.f90, fieldline_integrals.f90, fieldline_integrands.f90,
│   fieldline_labels.f90, make_fieldline.f90, surface_average.f90
├── shaing_callen_lib (src/shaing_callen/)
│   shaing_callen.f90, shaing_callen_integration.f90, shaing_callen_wrappers.f90
├── netcdf_lib (src/netcdf/)
│   netcdf.f90 - NetCDF output via type netcdf_t (backed by fortio)
├── precession_lib (src/precession/)
│   grid.f90 - Integration grid type and bounce/spline routines
│   precession.f90 - Magnetic well and turning-point helpers
│   saved_bounce_integrals.f90 - Per-fieldline bounce integrals of
│                        electric_rabe saved to file (type bounce_integrals_header_t)
├── error_handling_lib (src/error_handling/)
│   error_handling.f90 - Sanity-check helpers: failed_sanity_check(),
│                        reset_failed_check_counter(), did_fail_any_sanity_check(),
│                        set_unsafe_mode()
└── Top-level modules:
    coefficients.f90, deviation.f90, fit_functions.f90, read_file.f90
```

External dependencies (all fetched at build time): fortio, libneo, quadpack.
No system NetCDF library is required.

### Application flow (`app/main.f90`)

1. Read namelist configuration from `rabe.in`
2. Select the `boozer_field_t` loader from `field_type` (`vmec_nc` ->
   `init_from_vmec`, `booz_xform` -> `init_from_boozmn`, `chartmap` ->
   `init_from_chartmap`) and initialize from the `.nc` file named by `field_file`
3. Per flux surface: fix field, compute field lines, deviation factors, surface averages, optionally Shaing-Callen trapped fraction; if `unsafe_mode = .true.` and a sanity check fails, NaN-fill that surface's outputs instead of halting
4. Write results to `rabe.nc` (NetCDF)

### Key derived types

- `field_t` - Abstract base for magnetic field (in `field_base`)
- `boozer_field_t` - Field in Boozer coordinates (in `boozer_field`); load with `init_from_vmec(vmec_file, ...)`, `init_from_boozmn(boozmn_file)`, or `init_from_chartmap(chartmap_file)`
- `fieldline_t` - Field line representation (in `fieldline_mod`)
- `surface_average_t` - Surface-averaged quantities (in `surface_average_mod`)
- `netcdf_t` - NetCDF output handle (in `netcdf_mod`)

### Tests (`test/`)

- `test/unit/` - Unit tests (label: `quick`)
- `test/integration/` - Integration tests (label: `quick` or `slow`)
- `test/integration/vmec/` - VMEC field comparison tests
- `test/plot/` - Visualization tests; see `test/plot/README.md` for API and patterns
- `test/golden/` - Golden record regression tests of `rabe.x` and `electric_rabe` (label: `golden`)
- `test/restart/` - End-to-end test of resuming `electric_rabe` from its saved bounce
  integrals (label: `slow`)
- `test/helpers/` - Shared test utilities (analytical fields, mock fields, readers, plot helpers)
- `test/internal/` - Private tests with legacy NEO-2 comparisons (git submodule, not in public repo)

## Key Development Patterns

- **Naming**: Derived types use `typename_t` (e.g., `fieldline_t`, `boozer_field_t`)
- **Polymorphism**: Field calculations use abstract `field_t` base type. Fortran arrays cannot be polymorphic — use wrapper types with abstract members.
- **Module organization**: One focused responsibility per module
- **Test labels**: `quick`, `slow`, `plot`, `current`, `internal`, `golden`

## Input/Output

- **Input**: Namelist file `rabe.in` (namelist `rabe_config`)
- **Field data**: VMEC `.nc` files
- **Output**: `rabe.nc` (NetCDF)
- **`electric_rabe`**: reads `rabe.in` (namelists `rabe_config` and
  `precession_config`, single surface only), writes `electric_rabe.nc` and the
  bounce integrals to `saved_bounce_integrals.dat`. That file is plain text: a
  `# key = value` header with the inputs the bounce integrals depend on, then one
  fixed-width line per fieldline and grid point. Each fieldline is appended as soon
  as it is computed; on start, an existing file is checked against the input
  (mismatch stops the run), complete fieldlines are reused and only missing ones
  are computed. Only a clean build of the same commit reuses the file.
  The file is kept after the run as a cache. `n_eta` in
  `precession_config` sets the eta grid points of the bounce integrals (default 200).

## Code Style

- 88-character line limit
- 4-space indentation
- Out-of-source builds (use `build/` directory)

## Keeping Documentation in Sync

When making code changes, update the relevant documentation:
- **This file (`CLAUDE.md`)** — if adding/removing/renaming modules, libraries, derived types, make targets, test labels, or changing I/O formats
- **Subdirectory READMEs** — check for and update any README.md in directories you modify
