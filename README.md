# FT-IR-Species-Analysis-MATLAB-Code
Translated FORTRAN complier to MATLAB syntax for FT-IR Species Analysis

# FT-IR Species MATLAB Analysis Code — User Guide

---

## What You Need Before Running

### Required Files

| File Type | Description | Typical Name |
|-----------|-------------|--------------|
| FTIR data files (Input into the GUI is a folder. Path: "C:\FTIR_DATA\FT-IR Data") | Measured transmittance spectra | `*.1`, `*.txt`, or `*.dat` |
| HITRAN line files (Input into the GUI is a folder. Path: "C:\FTIR_DATA\HITRAN") | Reduced spectroscopic parameters | Example files: `NO.RED`, `N2O_2.RED`, |
| Partition function list (Input into the GUI is a text file. Path: "C:\FTIR_DATA\QUANT\Partition-Sums\Input_file_data.txt") | Points to PF data files | `Input_file_data.txt` |

***If Hitran files or Partition Sum folders & files are needed, please reach out to the repository owner for copies.***

### Required Software
- MATLAB R2016b or later
- No additional toolboxes required

### All .m Files Must Be On The MATLAB Path
Place all `.m` files in one folder and add it to the MATLAB path, or
run from that folder directly.

---

## Step-By-Step First Run

### 1. Open MATLAB and run:
```matlab
main_driver
```

### 2. Fill In The Parameter Fields

| Field | Description | Example Value |
|-------|-------------|---------------|
| Path Length (cm) | Physical length of gas cell | Example input: `5` |
| Pressure (atm) | Total gas pressure | Example input: `1.0` |
| Temperature (K) | Gas temperature in Kelvin | Example input: `296.2` |
| Min Wavenumber | Lower end of fitting range | Example input: `1850` |
| Max Wavenumber | Upper end of fitting range | Example input: `2206` |
| ILS Half-Width | Half-width of instrument sinc² function | Example input: `1.28565` |
| Fine Grid Step | Internal OD grid spacing | Example input: `0.015` |
| Initial Baseline a1 | Starting guess for baseline amplitude | Example input: `1.0` |
| Initial Baseline b1 | Starting guess for baseline slope | Example input: `-1e-7` |
| Summary Output File | Name of the results summary file | Example input: `results.O` |

### 3. Browse For FTIR Data Folder
Click Browse next to FTIR Data Folder and select the folder containing
your `.1` scan files. All files in this folder will be processed in sequence.

### 4. Browse For HITRAN Folder
Click Browse next to HITRAN Folder and select the folder containing your
reduced HITRAN files (`.RED` or `.dat`). The program will:
- Auto-detect all files in the folder (up to 6 molecules)
- Auto-count the number of lines in each file
- Populate the molecule table automatically

### 5. Fill In The Molecule Table
After browsing the HITRAN folder, complete each row:

| Column | What To Enter | Example |
|--------|---------------|---------|
| HITRAN File | Auto-filled — do not change | Example input: `NO.RED` |
| Mol ID | HITRAN molecule number | Example input: `8` for NO, `4` for N2O |
| Name | Short label for output files | Example input: `NO`, `N2O` |
| Init PP (atm) | Starting partial pressure guess | Example input: `0.1` |
| Num Lines | Auto-filled — do not change | Example input: `2964` |

**Common Molecule IDs (Refer to molecule_index_Version20.csv for molecule IDs 
   Path: "C:\FTIR_DATA\molecule_index_Version20.csv"):**

| ID | Molecule | ID | Molecule |
|----|----------|----|----------|
| 1 | H2O | 6 | CH4 |
| 2 | CO2 | 7 | O2 |
| 3 | O3 | 8 | NO |
| 4 | N2O | 9 | SO2 |
| 5 | CO | 10 | NO2 |

### 6. Browse For Partition Function File
Click Browse next to Partition Sum List File and select your
`Input_file_data.txt` file.

### 7. Browse For Output Folder
Click Browse next to Output Folder. All output files will be written here.

### 8. Click Proceed
The program closes the GUI and begins processing.

---

## Understanding The Output (Below is an example output)

### Command Window Output
```
IMX = 27733  (DVI=0.0150, DWV=1.28565)    <- confirm IMX is large
Loading partition functions...
QTofi: loaded 93 species x 705 temperature points.
Found 3 FTIR data file(s) to process.     <- number of scan files found

Processing scan 1 of 3: C:\...\scan001.1
  NDP = 277 data points in spectral range.
Initial chiSq = 2.1988725e+01             <- starting fit quality

Iter  NC  chiSq_before    chiSq_after     lambdaDamp
   1   0   2.199e+01       1.549e-01      1e-04   <- large drop = good
   2   0   1.549e-01       3.861e-04      1e-05   <- still improving
   3   0   3.861e-04       4.194e-05      1e-06
   4   0   4.194e-05       4.183e-05      1e-07
   5   0   4.183e-05       4.183e-05      1e-08
Converged at iteration 5
Final fitParams = 7.89e-05  9.60e-04  1.007e+00  3.18e-07
```

### Output Files Per Scan

**`scanname_coeff.txt`** — Fitted results:
```
CHI^2  =   0.4183161E-04
NO     =  0.960431E-03      <- partial pressure in atm
N2O    =  0.788652E-04
A      =  0.100727E+01      <- baseline amplitude
B      =  0.318142E-06      <- baseline slope
SIG(NO ) =  0.139703E-02    <- 1-sigma uncertainty
...
NO OF ITERATION IS = 5
```

**`scanname_transmittance.txt`** — Point-by-point comparison:
```
wavenumber   measured_T   model_T    residual
1850.0425    0.9931600    0.9920786  0.108141E-02
1851.3281    0.9951900    0.9956293 -0.439299E-03
...
```

**`results.O`** — One-line summary per scan:
```
  1     NO    0.9604E-01     N2O   0.7887E-02
  2     NO    0.9418E-01     N2O   0.7923E-02
  3     NO    0.9631E-01     N2O   0.7801E-02
```
Note: values are multiplied by 100 (percentage units).
Converged scans have 1 leading space.
Non-converged scans have 5 leading spaces — easy to identify.

---

## Batch Processing Multiple Scans

The program automatically handles any number of scan files:

1. Place all `.1` files in one folder
2. Browse to that folder in the GUI
3. Click Proceed
4. The program loops through every file and writes separate output files
   for each scan, plus one summary line per scan in `results.O`

There is no limit on the number of scan files.

---

## Choosing The Right Wavenumber Range

The fitting range should:
- Cover the absorption bands of all molecules being fitted
- Not be so wide that unrelated features confuse the baseline fit
- Include some baseline on both sides of each absorption band

**Common ranges:**

| Molecules | Suggested waveMin | Suggested waveMax |
|-----------|----------------|----------------|
| NO + N2O | 1850 | 2206 |
| CO₂ | 2250 | 2400 |
| CO | 2050 | 2200 |
| H2O | 1300 | 1900 |
| CH4 | 2900 | 3100 |

---

## Choosing Initial Partial Pressure Guesses

The initial guess affects how quickly the LM fitter converges:

| Situation | Recommended Init PP |
|-----------|---------------------|
| Unknown concentration | `0.1` atm (FORTRAN default) |
| Near-ambient trace gas | `1e-6` atm (sub-ppm starting point) |
| High concentration cell | `1e-3` to `1e-2` atm |
| Previously measured | Use last known value |

If the model shows deep absorption at the initial guess (T < 0.9 everywhere),
the initial PP is too large — reduce by 10× and try again.

If the model shows no absorption (T ≈ 1.0 everywhere), the initial PP is
too small or the molecule file/ID is wrong.

---

## Saving and Reusing Settings

Click Save Settings before Proceed to save all current GUI values to a
`.mat` file. To reuse settings in a future run, load the `.mat` file in
MATLAB before calling `main_driver`:

```matlab
load('my_settings.mat');
% Then run main_driver and manually re-enter values
```

---

## Tips For Best Results

1. **Temperature is in Kelvin in the GUI.**
   If your lab records temperature in Celsius, add 273.15 before entering.

2. **Molecule order in the table does not matter** — the program fits all
   molecules simultaneously regardless of order.

3. **If a molecule converges to 1e-10** (the positivity floor), it means
   either the molecule is not present in the sample or its absorption band
   is not covered by the chosen wavenumber range.
