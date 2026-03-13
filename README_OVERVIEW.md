# FTIR MATLAB Code — High-Level Overview

## What This Program Does

This FTIR MATLAB Code is a nonlinear least-squares spectral fitting program that retrieves
molecular concentrations from Fourier Transform Infrared (FTIR) transmittance
spectra. It compares a measured transmittance spectrum against a physics-based
model built from HITRAN spectroscopic line parameters, then adjusts molecular
partial pressures and a baseline function until the model matches the
measurement as closely as possible.

---

## How It Works — Step By Step

### Step 1 — User Input (GUI)
The program opens a graphical interface where the user provides:
- Sample conditions: path length, pressure, temperature
- Spectral range: min/max wavenumber, grid parameters fineGridStep and ilsHalfWidth
- File locations: FTIR data folder, HITRAN folder, partition function file,
  output folder
- Molecule table: one row per molecule with HITRAN file, molecule ID, name,
  and initial partial pressure guess

### Step 2 — Partition Function Loading (QTofi)
Reads a text file listing partition function data files for up to 93 molecular
species across 705 temperature points. These are used to correct HITRAN line
strengths from the reference temperature (296 K) to the actual sample
temperature.

### Step 3 — Optical Depth Array Construction (INPUT + OPTD)
For each molecule:
- Reads the reduced HITRAN line parameter file (.RED or .dat)
- For each spectral line within the fitting window (waveMin±45 cm⁻¹):
  - Applies temperature correction to line strength using partition functions
  - Computes pressure-broadened Lorentzian half-width
  - Calls OPTD to add each line's Lorentzian contribution to the fine optDepth grid
- The result is a 2D array optDepth(nGridPoints+1, nMolecules) where:
  - nGridPoints = fix((waveMax - waveMin + 60) / fineGridStep) grid points on the fine grid
  - fineGridStep = 0.015 cm⁻¹ typical fine grid spacing
  - nMolecules = number of molecules

### Step 4 — Forward Model (CONV + APND + TAUA)
For each measured data point at wavenumber V1:
- CONV integrates exp(-sum(optDepth*PP)) * baseline * ILS over ±30 cm⁻¹ window
- APND computes the sinc² Instrument Line Shape (ILS) at each offset
- TAUA evaluates the transmittance and its derivatives with respect to all
  fitted parameters
- The result is modeled transmittance T(V1) and derivatives dT/dfitParams(k)

### Step 5 — Nonlinear Fitting (MRQMIN + MRQCOF + GAUSSJ)
The Levenberg-Marquardt (LM) algorithm iteratively minimizes:

    CHI² = sum over all data points of (T_measured - T_model)²

- MRQMIN controls the LM step: accepts or rejects trial solutions,
  adjusts the damping parameter lambdaDamp
- MRQCOF builds the Hessian approximation (hessian matrix) and gradient
  (gradient vector) by accumulating Jacobian contributions over all data points
- GAUSSJ solves the linear system hessian * paramStep = gradient for the parameter step paramStep
- Convergence when |DCHISQ| ≤ 1e-7 or |DCHISQ|/CHISQ ≤ 1e-8

### Step 6 — Output (OUTPT)
For each scan, two output files are written:
- `*_coeff.txt` — fitted parameters, standard deviations, covariance matrix,
  iteration count
- `*_transmittance.txt` — wavenumber, measured T, model T, residual columns
A summary line is appended to the `results.O` summary file.

---

## Fitted Parameters

| Index | Parameter | Description |
|-------|-----------|-------------|
| fitParams(1..nMolecules) | Partial pressures (atm) | One per molecule |
| fitParams(nMolecules+1) | a1 (baseline amplitude) | Scales overall transmittance |
| fitParams(nMolecules+2) | b1 (baseline slope) | Linear wavenumber-dependent term |

Model transmittance at wavenumber V:

    T(V) = (a1 + b1*V) * integral[ exp(-sum_k(OD_k * PP_k)) * ILS(V-V') dV' ]

---

## Two Critical Grid Parameters

| Parameter | Example Value | Role |
|-----------|-----------|------|
| Fine Grid Step | 0.015 cm⁻¹ | Spacing of optDepth array — must resolve individual spectral lines |
| ILS Half-Width | 1.28565 cm⁻¹ | Half-width of sinc² instrument function — matches data spacing |

These are NOT interchangeable. fineGridStep controls the resolution of the internal
physics model. ilsHalfWidth controls the instrument broadening applied to the model
before comparing to measured data.

---

## File Structure

```
main_driver.m      — GUI, scan loop, LM iteration control
INPUT.m            — Read HITRAN files, build optDepth array
OPTD.m             — Add individual line Lorentzian to optDepth grid
INFR.m             — Convert wavenumber to grid index
RLFR.m             — Convert grid index to wavenumber
CONV.m             — Convolve optDepth with ILS, compute T and derivatives
APND.m             — Sinc² instrument line shape function
TAUA.m             — Evaluate transmittance and Jacobian at one point
MRQMIN.m           — Levenberg-Marquardt controller
MRQCOF.m           — Build Hessian and gradient matrices
GAUSSJ.m           — Gauss-Jordan linear system solver
INDAT.m            — Read measured FTIR transmittance data file
OUTPT.m            — Write coefficient and transmittance output files
QTofi.m            — Load partition function data
QofT_matlab.m      — Interpolate partition function at temperature T
```

---

## Batch Processing

Place multiple FTIR scan files (.1, .txt, or .dat) in the FTIR data folder.
The program automatically processes all files in sequence and writes:
- One `*_coeff.txt` per scan
- One `*_transmittance.txt` per scan
- One line per scan in the summary `results.O` file

Converged scans are written with a leading space.
Non-converged scans are written with leading spaces and a wider indent
so they can be identified in the summary file.

---

## Hardcoded Constants

| Constant | Value | Meaning |
|----------|-------|---------|
| PI1 | 3.14159 | Pi |
| XNL | 2.479×10¹⁹ | Loschmidt number (molecules/cm³/atm) |
| C2 | 1.438786 cm·K | Second radiation constant hc/k |
| Tref | 296.0 K | HITRAN reference temperature |
