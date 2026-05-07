# Graph Audit Analysis - COVID-19 Project

## Executive Summary
**Status: CRITICAL GAPS IDENTIFIED**

The `.tex` file references **18 graphs** but **only 1 exists** (the logo). The graphs need to be generated from optimal parameters extracted from the reference PDF.

---

## 1. GRAPHS REFERENCED IN .TEX FILE vs FILES PRESENT

### ✓ Present in Folder (1 file)
- `iiit_kalyani_logo.png` ✓

### ✗ Missing from Folder (17 files)
Required by `COVID_19_Modeling_uploaded_fixed.tex`:

1. `BifurCordi.png` - Bifurcation coordination
2. `R_beta_coord.png` - R vs beta coordination  
3. `R_vs_alpha.png` - R vs alpha parameter
4. `R_vs_b2Coord.png` - R vs b2 coordination
5. `R_vs_rho2Coord.png` - R vs rho2 coordination
6. `alpha_vs_beta_onR.png` - 2D: alpha vs beta on R0
7. `rho1_vs_rho2onR.png` - 2D: rho1 vs rho2 on R0
8. `A_vs_betaonR.png` - 2D: A vs beta on R0
9. `A_vs_donR.png` - 2D: A vs d on R0
10. `A_vs_rho2onR.png` - 2D: A vs rho2 on R0
11. `A_vs_rho1onR.png` - 2D: A vs rho1 on R0
12. `A_vs_b1onR.png` - 2D: A vs b1 on R0
13. `A_vs_b2.png` - 2D: A vs b2 on R0
14. `A_vs_alphaonR.png` - 2D: A vs alpha on R0
15. `A_vs_eta.png` - 2D: A vs eta on R0
16. `A_vs_conR.png` - 2D: A vs c on R0
17. `PRCC_Sensitivity.png` - PRCC sensitivity analysis

---

## 2. WHAT'S IN THE REFERENCE PDF

**File:** `Covid 19_Optimal Control.pdf`  
**Contains:** 24 figures from an S-I1-I2-H-R-V model with optimal control analysis

### Figures Available in PDF
- Time evolution plots
- DFE (Disease-Free Equilibrium) stability analysis
- Endemic equilibrium plots
- Bifurcation diagrams
- Sensitivity analysis (Tornado plots, heatmaps)
- Phase portraits
- R0 dependency plots
- Vaccination effect comparisons
- Transmission variation effects
- Optimal control impact visualizations
- PRCC (Partial Rank Correlation Coefficient) plots
- Control strategies comparison

---

## 3. MODEL MISMATCH ANALYSIS

### Current .tex File Model
**Model Type:** SEQIR (Susceptible-Exposed-Quarantined-Infected-Recovered)  
**Parameters:** A, beta, rho1, rho2, d, p, M, b1, b2, alpha, sigma, c, eta, delta

### PDF Reference Model
**Model Type:** SI1-I2-HR-V (with vaccination and hospitalization)  
**Parameters:** Lambda, beta1, beta2, d, eta, q1, q2, kappa, mu1, mu2, phi1, phi2, alpha, psi

### Result
❌ **Models use DIFFERENT parameters** → Graphs from PDF cannot directly populate current .tex file

---

## 4. ACTION REQUIRED

### Option A: Generate All Missing Graphs (Recommended)
Generate all 17 missing graphs using **optimal parameters from the PDF** for the SEQIR model:

1. **Extract optimal parameters** from `Covid 19_Optimal Control.pdf`
2. **Run sensitivity analysis** to generate:
   - 2D parameter plane graphs (A_vs_*, R_vs_*, alpha_vs_beta, etc.)
   - Bifurcation coordination plot
   - PRCC sensitivity analysis
3. **Verify all 18 graphs** are present before compiling .tex

### Option B: Update .tex File
Align the .tex file with the PDF's model (S-I1-I2-H-R-V) and use those 24 figures instead.

---

## 5. CURRENT PYTHON SCRIPTS STATUS

### Generate_Real_Graphs.py
- **Model:** SEQIR ✓
- **Output:** R0 calculations, trajectory plots
- **Status:** Missing the 2D parameter sensitivity plots needed for .tex

### Generate_Final_Figures.py
- **Model:** S-I1-I2-H-R-V (matches PDF)
- **Output:** Figures 2-8 (time evolution, stability, bifurcation, sensitivity, vaccination, transmission effects)
- **Status:** Could generate 7 figures from the PDF, but .tex wants different graphs

### covid_complete_plots.py
- **Status:** Not yet analyzed

---

## 6. RECOMMENDATIONS

✅ **Immediate Action:**
1. Identify the **optimal parameter set** from the PDF (usually in a table or results section)
2. Create a script that generates all 17 missing graphs using these optimal parameters
3. Verify parameter names match the SEQIR model implementation
4. Run the generation script to produce all PNG files
5. Compile the .tex file and verify all graphs render correctly

✅ **Documentation:**
- Add the optimal parameter values to `PARAMS` dictionary
- Add comments showing which PDF figure corresponds to which generated graph
- Version control the generation script with the parameter set used

---

## File Locations
- **Reference PDF:** `Covid 19_Optimal Control.pdf`
- **.tex file:** `COVID_19_Modeling_uploaded_fixed.tex`  
- **Python scripts:** `Generate_Real_Graphs.py`, `Generate_Final_Figures.py`
- **Missing image directory:** `c:\Users\mohda\Desktop\sem6\` (17 files needed)

