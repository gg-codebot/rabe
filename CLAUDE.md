# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

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
- `test/golden/` - Golden record regression test (label: `golden`)
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

## Code Style

- 88-character line limit
- 4-space indentation
- Out-of-source builds (use `build/` directory)

## Keeping Documentation in Sync

When making code changes, update the relevant documentation:
- **This file (`CLAUDE.md`)** — if adding/removing/renaming modules, libraries, derived types, make targets, test labels, or changing I/O formats
- **Subdirectory READMEs** — check for and update any README.md in directories you modify
