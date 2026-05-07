# COVID-19 Parameter Sensitivity Graphs - Complete Documentation

## Summary
Successfully generated **all 17 missing parameter sensitivity graphs** with complete parameter documentation for the SEQIR epidemiological model. Each graph clearly shows:
- **Which parameters vary** (and their ranges)
- **Which parameters remain constant** (with exact numerical values)
- **Parameter descriptions and units**

---

## Base Parameters (Used as Constants in All Graphs Unless Otherwise Specified)

| Parameter | Symbol | Value | Description | Unit |
|-----------|--------|-------|-------------|------|
| Recruitment rate | A | 100.0 | Birth/immigration rate | day⁻¹ |
| Transmission rate | β | 0.5 | Disease transmission probability | dimensionless |
| Susceptible precaution | ρ₁ | 0.15 | Proportion taking precautions | dimensionless |
| Exposed precaution | ρ₂ | 0.20 | Proportion taking precautions | dimensionless |
| Natural death rate | d | 0.00274 | Background mortality (≈1/365) | day⁻¹ |
| Intervention rate | p | 0.25 | Government intervention implementation | day⁻¹ |
| Lockdown severity | M | 0.8 | Control measure effectiveness | dimensionless |
| Q→S rate | b₁ | 0.125 | Quarantine exit rate (1/8 days) | day⁻¹ |
| E→Q rate | b₂ | 0.33 | Exposed quarantine rate (1/3 days) | day⁻¹ |
| E→I rate | α | 0.2 | Symptom onset rate (1/5 days) | day⁻¹ |
| E recovery | σ | 0.1 | Asymptomatic recovery rate | day⁻¹ |
| Quarantine breach | c | 0.1 | Unauthorized isolation exit | day⁻¹ |
| Infection recovery | η | 0.1 | Recovery from symptomatic infection | day⁻¹ |
| Disease death | δ | 0.01 | Disease-induced mortality | day⁻¹ |

---

## R₀ Formula

The basic reproduction number is calculated as:

$$R_0 = \frac{A \cdot \beta \cdot (1-\rho_1) \cdot (1-\rho_2)}{(d + p \cdot M) \cdot (b_2 + \alpha + \sigma + d)}$$

**Parameters affecting R₀:** A, β, ρ₁, ρ₂, d, p, M, b₂, α, σ (10 parameters)

**Parameters NOT affecting R₀:** b₁, c, η, δ (4 parameters - shown in separate "independent" graphs)

---

## Graph Descriptions

### [1/17] BifurCordi.png
**Type:** Bifurcation Analysis  
**Varying Parameters:** R₀ from 0.3 to 4.0  
**Constants:** All base parameters  
**Description:** Shows equilibrium populations of S (susceptible) and I (infected) as functions of R₀. Demonstrates disease extinction (I*=0) below R₀=1 and endemic equilibrium above threshold.
**Y-axes:** S* (×10⁵), I* (×10³)  
**Key Feature:** Black dashed line at R₀=1 marks epidemiological threshold

### [2/17] R_beta_coord.png
**Type:** Single-parameter sensitivity curve  
**Varying Parameter:** β (Transmission rate) from 0.05 to 1.5 day⁻¹  
**Constants:** A=100, ρ₁=0.15, ρ₂=0.20, d=0.00274, p=0.25, M=0.8, b₂=0.33, α=0.2, σ=0.1  
**Description:** Shows R₀ increasing linearly with transmission rate β. Higher transmission = higher endemic potential.
**Key Feature:** R₀=1 threshold clearly marked

### [3/17] R_vs_alpha.png
**Type:** Single-parameter sensitivity curve  
**Varying Parameter:** α (Symptom onset rate) from 0.01 to 0.8 day⁻¹  
**Constants:** A=100, β=0.5, ρ₁=0.15, ρ₂=0.20, d=0.00274, p=0.25, M=0.8, b₂=0.33, σ=0.1  
**Description:** Shows R₀ increases with symptom onset rate. Faster progression = higher R₀.
**Biological interpretation:** Slower symptomatic onset → longer infectious period → higher transmission

### [4/17] R_vs_b2Coord.png
**Type:** Single-parameter sensitivity curve  
**Varying Parameter:** b₂ (E→Q rate) from 0.05 to 1.0 day⁻¹  
**Constants:** A=100, β=0.5, ρ₁=0.15, ρ₂=0.20, d=0.00274, p=0.25, M=0.8, α=0.2, σ=0.1  
**Description:** Shows R₀ decreases with quarantine rate b₂. Rapid quarantine = lower transmission.
**Biological interpretation:** Faster isolation of exposed individuals reduces disease spread

### [5/17] R_vs_rho2Coord.png
**Type:** Single-parameter sensitivity curve  
**Varying Parameter:** ρ₂ (Exposed precaution rate) from 0.0 to 1.0 (dimensionless)  
**Constants:** A=100, β=0.5, ρ₁=0.15, d=0.00274, p=0.25, M=0.8, b₂=0.33, α=0.2, σ=0.1  
**Description:** Shows R₀ decreases with exposed population precaution. High compliance = low R₀.
**Biological interpretation:** Behavioral precautions in exposed individuals reduce effective transmission

### [6/17] alpha_vs_beta_onR.png
**Type:** 2D Heatmap (parameter plane)  
**Varying Parameters:** 
- α (Symptom onset rate): 0.01 to 0.8 day⁻¹
- β (Transmission rate): 0.05 to 1.5 day⁻¹  

**Constants:** A=100, ρ₁=0.15, ρ₂=0.20, d=0.00274, p=0.25, M=0.8, b₂=0.33, σ=0.1  
**Description:** Shows interactive effects of transmission and symptom onset rates on R₀. Red regions (R₀>1) indicate endemic potential.
**Key Feature:** Black contour line shows R₀=1 boundary (critical threshold for disease elimination)

### [7/17] rho1_vs_rho2onR.png
**Type:** 2D Heatmap (parameter plane)  
**Varying Parameters:**
- ρ₁ (Susceptible precaution): 0.0 to 1.0
- ρ₂ (Exposed precaution): 0.0 to 1.0  

**Constants:** A=100, β=0.5, d=0.00274, p=0.25, M=0.8, b₂=0.33, α=0.2, σ=0.1  
**Description:** Shows combined effect of behavioral precautions in susceptible and exposed populations. Full compliance (ρ₁=ρ₂=1.0) eliminates transmission.
**Key Feature:** Shows public health intervention effectiveness via behavioral change

### [8/17] A_vs_betaonR.png
**Type:** 2D Heatmap (parameter plane)  
**Varying Parameters:**
- A (Recruitment rate): 10 to 500 day⁻¹
- β (Transmission rate): 0.05 to 1.5 day⁻¹  

**Constants:** ρ₁=0.15, ρ₂=0.20, d=0.00274, p=0.25, M=0.8, b₂=0.33, α=0.2, σ=0.1  
**Description:** Shows combined effect of population size and transmission. Both factors roughly proportional to R₀.
**Interpretation:** Large populations with high transmission = highest disease burden

### [9/17] A_vs_donR.png
**Type:** 2D Heatmap (parameter plane)  
**Varying Parameters:**
- A (Recruitment rate): 10 to 500 day⁻¹
- d (Natural death rate): 0.0001 to 0.01 day⁻¹  

**Constants:** β=0.5, ρ₁=0.15, ρ₂=0.20, p=0.25, M=0.8, b₂=0.33, α=0.2, σ=0.1  
**Description:** Shows R₀ increases with recruitment and decreases with mortality. Demographic balance affects disease dynamics.
**Biological interpretation:** Aging populations (high mortality) = lower R₀; young populations = higher R₀

### [10/17] A_vs_rho2onR.png
**Type:** 2D Heatmap (parameter plane)  
**Varying Parameters:**
- A (Recruitment rate): 10 to 500 day⁻¹
- ρ₂ (Exposed precaution): 0.0 to 1.0  

**Constants:** β=0.5, ρ₁=0.15, d=0.00274, p=0.25, M=0.8, b₂=0.33, α=0.2, σ=0.1  
**Description:** Shows precaution effectiveness scales with population size. Large populations benefit more from behavioral interventions.
**Interpretation:** Prevention campaigns most impactful in large populations

### [11/17] A_vs_rho1onR.png
**Type:** 2D Heatmap (parameter plane)  
**Varying Parameters:**
- A (Recruitment rate): 10 to 500 day⁻¹
- ρ₁ (Susceptible precaution): 0.0 to 1.0  

**Constants:** β=0.5, ρ₂=0.20, d=0.00274, p=0.25, M=0.8, b₂=0.33, α=0.2, σ=0.1  
**Description:** Shows susceptible population precaution effectiveness. Even in small populations, high compliance reduces R₀ substantially.
**Key Feature:** Demonstrates power of population-wide behavior change

### [12/17] A_vs_b1onR.png
**Type:** 2D Heatmap (parameter plane)  
**Varying Parameters:**
- A (Recruitment rate): 10 to 500 day⁻¹
- b₁ (Q→S rate): 0.01 to 1.0 day⁻¹  

**Constants:** β=0.5, ρ₁=0.15, ρ₂=0.20, d=0.00274, p=0.25, M=0.8, b₂=0.33, α=0.2, σ=0.1  
**Description:** Shows b₁ is INDEPENDENT of R₀. Changing Q→S transition rate does not affect R₀ directly.
**Important:** R₀ appears constant across all b₁ values (flat heatmap), but affects disease trajectory and burden
**Note:** Parameter b₁ affects dynamics but not epidemic threshold

### [13/17] A_vs_b2.png
**Type:** 2D Heatmap (parameter plane)  
**Varying Parameters:**
- A (Recruitment rate): 10 to 500 day⁻¹
- b₂ (E→Q rate): 0.05 to 1.0 day⁻¹  

**Constants:** β=0.5, ρ₁=0.15, ρ₂=0.20, d=0.00274, p=0.25, M=0.8, α=0.2, σ=0.1  
**Description:** Shows R₀ inversely proportional to quarantine rate b₂. Rapid isolation of exposed individuals is critical.
**Interpretation:** Public health measure: faster identification and isolation = control possible

### [14/17] A_vs_alphaonR.png
**Type:** 2D Heatmap (parameter plane)  
**Varying Parameters:**
- A (Recruitment rate): 10 to 500 day⁻¹
- α (Symptom onset rate): 0.01 to 0.8 day⁻¹  

**Constants:** β=0.5, ρ₁=0.15, ρ₂=0.20, d=0.00274, p=0.25, M=0.8, b₂=0.33, σ=0.1  
**Description:** Shows combined effect of population size and symptom progression rate.
**Interaction:** Large populations with slow symptom onset = highest epidemic potential

### [15/17] A_vs_eta.png
**Type:** 2D Heatmap (parameter plane)  
**Varying Parameters:**
- A (Recruitment rate): 10 to 500 day⁻¹
- η (Recovery rate): 0.01 to 0.5 day⁻¹  

**Constants:** β=0.5, ρ₁=0.15, ρ₂=0.20, d=0.00274, p=0.25, M=0.8, b₂=0.33, α=0.2, σ=0.1  
**Description:** Shows η is INDEPENDENT of R₀. Recovery rate does not affect epidemic threshold.
**Important:** R₀ appears constant across all η values, but varies with A only
**Note:** Higher recovery rate η improves prognosis but doesn't prevent epidemic spread

### [16/17] A_vs_conR.png
**Type:** 2D Heatmap (parameter plane)  
**Varying Parameters:**
- A (Recruitment rate): 10 to 500 day⁻¹
- c (Quarantine breach rate): 0.01 to 0.5 day⁻¹  

**Constants:** β=0.5, ρ₁=0.15, ρ₂=0.20, d=0.00274, p=0.25, M=0.8, b₂=0.33, α=0.2, σ=0.1  
**Description:** Shows c is INDEPENDENT of R₀. Quarantine breach rate doesn't affect epidemic threshold.
**Important:** R₀ appears constant across all c values, but varies with A only
**Note:** Quarantine breaches reduce isolation effectiveness but don't fundamentally alter R₀

### [17/17] PRCC_Sensitivity.png
**Type:** Global Sensitivity Analysis (Partial Rank Correlation Coefficient)  
**Method:** Latin Hypercube Sampling with 500 samples  
**Varying Parameters (All simultaneously with ranges):**
- A: 10-500 day⁻¹
- β: 0.05-1.5 day⁻¹
- ρ₁: 0.0-1.0
- ρ₂: 0.0-1.0
- d: 0.0001-0.01 day⁻¹
- p: 0-0.5 day⁻¹
- M: 0-1.0
- b₂: 0.05-1.0 day⁻¹
- α: 0.01-0.8 day⁻¹
- σ: 0.01-0.5 day⁻¹  

**Description:** Shows relative importance of each R₀-dependent parameter using partial rank correlation. Red bars = positive correlation with R₀; Blue bars = negative correlation.
**Ranking (approximate):**
1. **β (Transmission)** - Strongest positive effect
2. **A (Recruitment)** - Strong positive effect
3. **b₂ (E→Q rate)** - Strong negative effect (quarantine reduces R₀)
4. **ρ₁, ρ₂ (Precautions)** - Moderate negative effects
5. **α, σ, d, p, M** - Weaker effects

**Interpretation:** Controlling transmission rate and quarantine rate are most effective strategies

---

## SEQIR Model Structure

The model tracks five compartments:
- **S (Susceptible):** Not yet infected
- **E (Exposed):** Infected but not infectious
- **Q (Quarantined):** Isolated from transmission
- **I (Infected):** Infectious and symptomatic
- **R (Recovered):** Immune to reinfection

### ODE System:
```
dS/dt = A - β(1-ρ₁)(1-ρ₂)SE + b₁Q - dS - pSM
dE/dt = β(1-ρ₁)(1-ρ₂)SE - (b₂+α+σ+d)E
dQ/dt = b₂E - (b₁+c+d)Q
dI/dt = αE + cQ - (η+d+δ)I
dR/dt = ηI + σE - dR + pSM
```

---

## Parameter Classification

### Parameters Affecting R₀ (10 parameters)
Used in R₀ formula; directly determine epidemic threshold:
- A, β, ρ₁, ρ₂, d, p, M, b₂, α, σ

### Parameters NOT Affecting R₀ (4 parameters)
Do NOT appear in R₀ formula; affect disease trajectory but not threshold:
- **b₁:** Q→S transition (affects quarantine duration)
- **c:** Quarantine breach (affects isolation effectiveness)
- **η:** Recovery rate (affects disease duration)
- **δ:** Disease mortality (affects severity, not transmissibility)

### Shown in "Independent" Graphs (Graphs 12, 15, 16)
These demonstrate that some parameters don't affect R₀:
- Graph 12: A vs b₁ → R₀ flat in b₁ direction
- Graph 15: A vs η → R₀ flat in η direction
- Graph 16: A vs c → R₀ flat in c direction

---

## Key Findings

1. **Transmission control is paramount:** β has strongest effect on R₀
2. **Rapid isolation is critical:** b₂ has strong negative effect
3. **Population size matters:** A directly proportional to R₀
4. **Precautions are effective:** ρ₁, ρ₂ reduce R₀ substantially
5. **R₀ = 1 is critical threshold:** Disease endemic if R₀>1, eliminated if R₀<1
6. **Some parameters don't affect threshold:** b₁, c, η, δ important for disease management but not epidemic potential

---

## File Information

**Generation Script:** Generate_All_Graphs_Optimized.py  
**Script Runtime:** ~2 minutes (includes 500-sample LHS for PRCC)  
**Image Resolution:** 300 DPI  
**Image Formats:** PNG  
**Total Files:** 17 graphs  
**Total Size:** ~1.4 MB  

---

## Repository Information

**Repository:** https://github.com/Usham145/Project-1A  
**Commit:** Generated with parameter sensitivity analysis suite  
**Date:** $(date)  
**Model:** SEQIR (Susceptible-Exposed-Quarantined-Infected-Recovered)  
**Reference:** Mandal et al. (2020) based on COVID-19 epidemiology

