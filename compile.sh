# Constant-Returns-Cash-Model

Fortran implementation of a constant-returns-to-scale (CRS) dynamic investment model with cash holdings. The model is estimated via simulated method of moments (SMM).

## Model Overview

Firms use capital in a CRS technology to generate operating income. Each period they choose investment and cash holdings subject to adjustment costs and external financing costs. The model is solved via value function iteration and simulated to produce moments matched to Compustat data.

## Folder Structure

```
fortran/
├── Input/
│   └── estfil.txt           # Parameter values for evaluation
├── Outputcomp/              # Model output
│   ├── policies.txt         # Full policy function output
│   ├── c.txt, i.txt, d.txt  # Cash, investment, dividend policies
│   ├── v.txt                # Value function
│   ├── cg.txt, zg.txt       # State grids (cash, productivity)
│   ├── statespaces.txt      # Grid specifications
│   ├── simdata.txt          # Simulated panel data
│   ├── moments.txt          # Computed moments
│   └── *.png                # Policy function plots
└── *.f90                    # Source files (see below)
```

## Source Files

### Modules (`crs_mod.f90`)
Contains all module definitions:
- `datatype`: Double precision and data structure definitions
- `sizes`: Grid dimensions, simulation parameters, output controls
- `globals`: Global variables (simulated data, random draws)
- `pickmoments`: Moment selection and naming
- `myseed`: Random number seed initialization

### Core Model Files

| File | Description |
|------|-------------|
| `crs_solve.f90` | Value function iteration (Bellman equation solver) |
| `modelfunctions.f90` | Grid construction, Bellman setup/maximization, policy output |
| `simmodel.f90` | Panel simulation using solved policy functions |
| `makemoments.f90` | Computes moments from simulated panel |
| `shock_draw.f90` | Generates random draws for simulation |

### Driver Programs

| File | Description |
|------|-------------|
| `evaluate_crs.f90` | Evaluates model at given parameters, reports moments |
| `compstat.f90` | Comparative statics analysis |

## Compilation

```bash
./compile.sh          # Release build
./compile.sh debug    # Debug build with bounds checking
```

Requires `gfortran`. The executable is named `crs_model`.

## Usage

```bash

# Run the model
./crs_model
```

Parameters are read from `Input/estfil.txt` if it exists, otherwise defaults are used.

## Key Parameters

| Parameter | Symbol | Description |
|-----------|--------|-------------|
| `psi` | γ | Investment adjustment cost |
| `delta` | δ | Depreciation rate |
| `mu` | μ | Mean of log productivity |
| `rho` | ρ | Persistence of productivity |
| `sig_z` | σ | Std. dev. of productivity shocks |

## Configuration (`crs_mod.f90`)

Key settings in the `sizes` module:
- `nc = 81`: Cash state grid points
- `ncp = 201`: Cash policy grid points
- `ni = 201`: Investment grid points
- `nz = 21`: Productivity grid points
- `nFirms = 1000`: Number of simulated firms
- `nYears = 110`: Simulation length (100 burn-in + 10 sample)
- `printfiles = 0/1`: Enable output file writing
- `verbose = 0/1`: Enable iteration output

## Output Files

When `printfiles=1`:
- `policies.txt`: Human-readable policy functions
- `c.txt`, `i.txt`, `d.txt`, `v.txt`: Policy matrices (rows=productivity, cols=cash states)
- `simdata.txt`: Simulated panel (firm, year, value, cash, investment, dividends, productivity)
- `moments.txt`: Summary of parameters and computed moments

