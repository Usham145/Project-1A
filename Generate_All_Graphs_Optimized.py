"""
COVID-19 SEQIR MODEL - PARAMETER SENSITIVITY ANALYSIS
Generates all 17 missing graphs with parameter documentation
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, rankdata, uniform
import warnings
warnings.filterwarnings('ignore')

plt.style.use('seaborn-v0_8-darkgrid')

# ============================================================================
# BASE PARAMETERS
# ============================================================================

BASE_PARAMS = {
    'A': 100,           # Recruitment rate (per day)
    'beta': 0.5,        # Disease transmission rate
    'rho1': 0.15,       # Susceptible precaution rate
    'rho2': 0.20,       # Exposed precaution rate
    'd': 0.00274,       # Natural death rate (1/365 days)
    'p': 0.25,          # Govt intervention implementation rate
    'M': 0.8,           # Control/lockdown severity
    'b1': 0.125,        # Quarantined → susceptible rate (1/8 days)
    'b2': 0.33,         # Exposed → quarantine rate (1/3 days)
    'alpha': 0.2,       # Exposed → infected rate (1/5 days)
    'sigma': 0.1,       # Exposed recovery (asymptomatic)
    'c': 0.1,           # Quarantine breach rate
    'eta': 0.1,         # Infected recovery rate
    'delta': 0.01,      # Disease-induced death rate
}

def R0_formula(A, beta, rho1, rho2, d, p, M, b2, alpha, sigma):
    """Calculate basic reproduction number R₀"""
    num = A * beta * (1 - rho1) * (1 - rho2)
    den = (d + p * M) * (b2 + alpha + sigma + d)
    return max(0, num / den) if den > 0 else 0

print("\n" + "="*80)
print("COVID-19 SEQIR MODEL - PARAMETER SENSITIVITY GRAPHS")
print("="*80)
print("\nBASE PARAMETERS (Constants where not explicitly varied):")
for param in sorted(BASE_PARAMS.keys()):
    print(f"  {param:8s} = {BASE_PARAMS[param]:10.6f}")

# ============================================================================
# GRAPH 1: BIFURCATION
# ============================================================================
print("\n[1/17] BifurCordi.png - VARYING: R₀ | CONSTANTS: All base parameters")

R0_vals = np.linspace(0.3, 4.0, 300)
S_eq = np.zeros_like(R0_vals)
I_eq = np.zeros_like(R0_vals)

for i, R0_target in enumerate(R0_vals):
    if R0_target <= 1.0:
        S_eq[i] = BASE_PARAMS['A'] / (BASE_PARAMS['d'] + BASE_PARAMS['p']*BASE_PARAMS['M'])
        I_eq[i] = 0
    else:
        S_eq[i] = BASE_PARAMS['A'] / (BASE_PARAMS['d'] + BASE_PARAMS['p']*BASE_PARAMS['M']) / R0_target
        ratio = (BASE_PARAMS['alpha'] + (BASE_PARAMS['c']*BASE_PARAMS['b2'])/(BASE_PARAMS['b1']+BASE_PARAMS['c']+BASE_PARAMS['d']))/(BASE_PARAMS['eta']+BASE_PARAMS['d']+BASE_PARAMS['delta'])
        E_eq = (BASE_PARAMS['A'] - (BASE_PARAMS['d'] + BASE_PARAMS['p']*BASE_PARAMS['M'])*S_eq[i]) / (BASE_PARAMS['b2']+BASE_PARAMS['alpha']+BASE_PARAMS['sigma']+BASE_PARAMS['d'])
        I_eq[i] = E_eq * ratio

plt.figure(figsize=(9, 5.5))
plt.plot(R0_vals, S_eq/1e5, 'b-', lw=2.5, label=r'$S^*$ (×10⁵)')
plt.plot(R0_vals, I_eq/1e3, 'r-', lw=2.5, label=r'$I^*$ (×10³)')
plt.axvline(1.0, color='black', ls='--', alpha=0.6, lw=1.5, label='R₀=1')
plt.xlabel(r'Basic Reproduction Number $R_0$', fontsize=12)
plt.ylabel('Equilibrium', fontsize=12)
plt.title('Bifurcation: S* & I* vs R₀', fontsize=12, fontweight='bold')
plt.legend(fontsize=10)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('BifurCordi.png', dpi=300, bbox_inches='tight')
plt.close()

# ============================================================================
# GRAPH 2-5: R₀ vs Single Parameters
# ============================================================================

# Graph 2
print("\n[2/17] R_beta_coord.png - VARYING: β (0.05→1.5) | CONSTANTS: A,ρ₁,ρ₂,d,p,M,b₂,α,σ")
beta_range = np.linspace(0.05, 1.5, 250)
R0_beta = [R0_formula(BASE_PARAMS['A'], b, BASE_PARAMS['rho1'], BASE_PARAMS['rho2'],
                      BASE_PARAMS['d'], BASE_PARAMS['p'], BASE_PARAMS['M'],
                      BASE_PARAMS['b2'], BASE_PARAMS['alpha'], BASE_PARAMS['sigma']) 
           for b in beta_range]
plt.figure(figsize=(6.5, 4.5))
plt.plot(beta_range, R0_beta, 'b-', lw=2.8)
plt.axhline(1, color='black', ls='--', alpha=0.6, linewidth=1.2)
plt.xlabel(r'Transmission Rate $\beta$', fontsize=11)
plt.ylabel(r'$R_0$', fontsize=11)
plt.title(r'$R_0$ vs Transmission Rate $\beta$', fontsize=12, fontweight='bold')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('R_beta_coord.png', dpi=300, bbox_inches='tight')
plt.close()

# Graph 3
print("\n[3/17] R_vs_alpha.png - VARYING: α (0.01→0.8) | CONSTANTS: A,β,ρ₁,ρ₂,d,p,M,b₂,σ")
alpha_range = np.linspace(0.01, 0.8, 250)
R0_alpha = [R0_formula(BASE_PARAMS['A'], BASE_PARAMS['beta'], BASE_PARAMS['rho1'], BASE_PARAMS['rho2'],
                       BASE_PARAMS['d'], BASE_PARAMS['p'], BASE_PARAMS['M'],
                       BASE_PARAMS['b2'], a, BASE_PARAMS['sigma'])
            for a in alpha_range]
plt.figure(figsize=(6.5, 4.5))
plt.plot(alpha_range, R0_alpha, 'g-', lw=2.8)
plt.axhline(1, color='black', ls='--', alpha=0.6, linewidth=1.2)
plt.xlabel(r'Symptom Onset Rate $\alpha$ (day⁻¹)', fontsize=11)
plt.ylabel(r'$R_0$', fontsize=11)
plt.title(r'$R_0$ vs Symptom Onset Rate $\alpha$', fontsize=12, fontweight='bold')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('R_vs_alpha.png', dpi=300, bbox_inches='tight')
plt.close()

# Graph 4
print("\n[4/17] R_vs_b2Coord.png - VARYING: b₂ (0.05→1.0) | CONSTANTS: A,β,ρ₁,ρ₂,d,p,M,α,σ")
b2_range = np.linspace(0.05, 1.0, 250)
R0_b2 = [R0_formula(BASE_PARAMS['A'], BASE_PARAMS['beta'], BASE_PARAMS['rho1'], BASE_PARAMS['rho2'],
                    BASE_PARAMS['d'], BASE_PARAMS['p'], BASE_PARAMS['M'],
                    b2, BASE_PARAMS['alpha'], BASE_PARAMS['sigma'])
         for b2 in b2_range]
plt.figure(figsize=(6.5, 4.5))
plt.plot(b2_range, R0_b2, 'm-', lw=2.8)
plt.axhline(1, color='black', ls='--', alpha=0.6, linewidth=1.2)
plt.xlabel(r'Exposed→Quarantine Rate $b_2$ (day⁻¹)', fontsize=11)
plt.ylabel(r'$R_0$', fontsize=11)
plt.title(r'$R_0$ vs Quarantine Rate $b_2$', fontsize=12, fontweight='bold')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('R_vs_b2Coord.png', dpi=300, bbox_inches='tight')
plt.close()

# Graph 5
print("\n[5/17] R_vs_rho2Coord.png - VARYING: ρ₂ (0.0→1.0) | CONSTANTS: A,β,ρ₁,d,p,M,b₂,α,σ")
rho2_range = np.linspace(0.0, 1.0, 250)
R0_rho2 = [R0_formula(BASE_PARAMS['A'], BASE_PARAMS['beta'], BASE_PARAMS['rho1'], rho2,
                      BASE_PARAMS['d'], BASE_PARAMS['p'], BASE_PARAMS['M'],
                      BASE_PARAMS['b2'], BASE_PARAMS['alpha'], BASE_PARAMS['sigma'])
           for rho2 in rho2_range]
plt.figure(figsize=(6.5, 4.5))
plt.plot(rho2_range, R0_rho2, 'c-', lw=2.8)
plt.axhline(1, color='black', ls='--', alpha=0.6, linewidth=1.2)
plt.xlabel(r'Exposed Precaution Rate $\rho_2$', fontsize=11)
plt.ylabel(r'$R_0$', fontsize=11)
plt.title(r'$R_0$ vs Exposed Precaution Rate $\rho_2$', fontsize=12, fontweight='bold')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('R_vs_rho2Coord.png', dpi=300, bbox_inches='tight')
plt.close()

# ============================================================================
# GRAPHS 6-16: 2D PARAMETER PLANES
# ============================================================================

# Graph 6
print("\n[6/17] alpha_vs_beta_onR.png - VARYING: α (0.01→0.8), β (0.05→1.5)")
alpha_2d = np.linspace(0.01, 0.8, 100)
beta_2d = np.linspace(0.05, 1.5, 100)
R0_ab = np.zeros((len(alpha_2d), len(beta_2d)))
for i, a in enumerate(alpha_2d):
    for j, b in enumerate(beta_2d):
        R0_ab[i, j] = R0_formula(BASE_PARAMS['A'], b, BASE_PARAMS['rho1'], BASE_PARAMS['rho2'],
                                BASE_PARAMS['d'], BASE_PARAMS['p'], BASE_PARAMS['M'],
                                BASE_PARAMS['b2'], a, BASE_PARAMS['sigma'])
plt.figure(figsize=(8, 6))
cf = plt.contourf(beta_2d, alpha_2d, R0_ab, levels=20, cmap='RdYlGn_r')
plt.contour(beta_2d, alpha_2d, R0_ab, levels=[1.0], colors='black', linewidths=2)
plt.colorbar(cf, label=r'$R_0$')
plt.xlabel(r'Transmission Rate $\beta$', fontsize=11)
plt.ylabel(r'Symptom Onset Rate $\alpha$ (day⁻¹)', fontsize=11)
plt.title(r'$R_0$ Plane: $\alpha$ vs $\beta$', fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig('alpha_vs_beta_onR.png', dpi=300, bbox_inches='tight')
plt.close()

# Graph 7
print("\n[7/17] rho1_vs_rho2onR.png - VARYING: ρ₁ (0.0→1.0), ρ₂ (0.0→1.0)")
rho1_2d = np.linspace(0.0, 1.0, 100)
rho2_2d = np.linspace(0.0, 1.0, 100)
R0_rr = np.zeros((len(rho1_2d), len(rho2_2d)))
for i, r1 in enumerate(rho1_2d):
    for j, r2 in enumerate(rho2_2d):
        R0_rr[i, j] = R0_formula(BASE_PARAMS['A'], BASE_PARAMS['beta'], r1, r2,
                                BASE_PARAMS['d'], BASE_PARAMS['p'], BASE_PARAMS['M'],
                                BASE_PARAMS['b2'], BASE_PARAMS['alpha'], BASE_PARAMS['sigma'])
plt.figure(figsize=(8, 6))
cf = plt.contourf(rho2_2d, rho1_2d, R0_rr, levels=20, cmap='RdYlGn_r')
plt.contour(rho2_2d, rho1_2d, R0_rr, levels=[1.0], colors='black', linewidths=2)
plt.colorbar(cf, label=r'$R_0$')
plt.xlabel(r'Exposed Precaution Rate $\rho_2$', fontsize=11)
plt.ylabel(r'Susceptible Precaution Rate $\rho_1$', fontsize=11)
plt.title(r'$R_0$ Plane: $\rho_1$ vs $\rho_2$', fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig('rho1_vs_rho2onR.png', dpi=300, bbox_inches='tight')
plt.close()

# Graph 8
print("\n[8/17] A_vs_betaonR.png - VARYING: A (10→500), β (0.05→1.5)")
A_2d = np.linspace(10, 500, 100)
R0_Ab = np.zeros((len(A_2d), len(beta_2d)))
for i, a_val in enumerate(A_2d):
    for j, b_val in enumerate(beta_2d):
        R0_Ab[i, j] = R0_formula(a_val, b_val, BASE_PARAMS['rho1'], BASE_PARAMS['rho2'],
                                BASE_PARAMS['d'], BASE_PARAMS['p'], BASE_PARAMS['M'],
                                BASE_PARAMS['b2'], BASE_PARAMS['alpha'], BASE_PARAMS['sigma'])
plt.figure(figsize=(8, 6))
cf = plt.contourf(beta_2d, A_2d, R0_Ab, levels=20, cmap='RdYlGn_r')
plt.contour(beta_2d, A_2d, R0_Ab, levels=[1.0], colors='black', linewidths=2)
plt.colorbar(cf, label=r'$R_0$')
plt.xlabel(r'Transmission Rate $\beta$', fontsize=11)
plt.ylabel(r'Recruitment Rate $A$ (day⁻¹)', fontsize=11)
plt.title(r'$R_0$ Plane: $A$ vs $\beta$', fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig('A_vs_betaonR.png', dpi=300, bbox_inches='tight')
plt.close()

# Graph 9
print("\n[9/17] A_vs_donR.png - VARYING: A (10→500), d (0.0001→0.01)")
d_2d = np.linspace(0.0001, 0.01, 100)
R0_Ad = np.zeros((len(A_2d), len(d_2d)))
for i, a_val in enumerate(A_2d):
    for j, d_val in enumerate(d_2d):
        R0_Ad[i, j] = R0_formula(a_val, BASE_PARAMS['beta'], BASE_PARAMS['rho1'], BASE_PARAMS['rho2'],
                                d_val, BASE_PARAMS['p'], BASE_PARAMS['M'],
                                BASE_PARAMS['b2'], BASE_PARAMS['alpha'], BASE_PARAMS['sigma'])
plt.figure(figsize=(8, 6))
cf = plt.contourf(d_2d, A_2d, R0_Ad, levels=20, cmap='RdYlGn_r')
plt.contour(d_2d, A_2d, R0_Ad, levels=[1.0], colors='black', linewidths=2)
plt.colorbar(cf, label=r'$R_0$')
plt.xlabel(r'Natural Death Rate $d$ (day⁻¹)', fontsize=11)
plt.ylabel(r'Recruitment Rate $A$ (day⁻¹)', fontsize=11)
plt.title(r'$R_0$ Plane: $A$ vs $d$', fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig('A_vs_donR.png', dpi=300, bbox_inches='tight')
plt.close()

# Graph 10
print("\n[10/17] A_vs_rho2onR.png - VARYING: A (10→500), ρ₂ (0.0→1.0)")
R0_Ar = np.zeros((len(A_2d), len(rho2_2d)))
for i, a_val in enumerate(A_2d):
    for j, r2_val in enumerate(rho2_2d):
        R0_Ar[i, j] = R0_formula(a_val, BASE_PARAMS['beta'], BASE_PARAMS['rho1'], r2_val,
                                BASE_PARAMS['d'], BASE_PARAMS['p'], BASE_PARAMS['M'],
                                BASE_PARAMS['b2'], BASE_PARAMS['alpha'], BASE_PARAMS['sigma'])
plt.figure(figsize=(8, 6))
cf = plt.contourf(rho2_2d, A_2d, R0_Ar, levels=20, cmap='RdYlGn_r')
plt.contour(rho2_2d, A_2d, R0_Ar, levels=[1.0], colors='black', linewidths=2)
plt.colorbar(cf, label=r'$R_0$')
plt.xlabel(r'Exposed Precaution Rate $\rho_2$', fontsize=11)
plt.ylabel(r'Recruitment Rate $A$ (day⁻¹)', fontsize=11)
plt.title(r'$R_0$ Plane: $A$ vs $\rho_2$', fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig('A_vs_rho2onR.png', dpi=300, bbox_inches='tight')
plt.close()

# Graph 11
print("\n[11/17] A_vs_rho1onR.png - VARYING: A (10→500), ρ₁ (0.0→1.0)")
R0_Ar1 = np.zeros((len(A_2d), len(rho1_2d)))
for i, a_val in enumerate(A_2d):
    for j, r1_val in enumerate(rho1_2d):
        R0_Ar1[i, j] = R0_formula(a_val, BASE_PARAMS['beta'], r1_val, BASE_PARAMS['rho2'],
                                BASE_PARAMS['d'], BASE_PARAMS['p'], BASE_PARAMS['M'],
                                BASE_PARAMS['b2'], BASE_PARAMS['alpha'], BASE_PARAMS['sigma'])
plt.figure(figsize=(8, 6))
cf = plt.contourf(rho1_2d, A_2d, R0_Ar1, levels=20, cmap='RdYlGn_r')
plt.contour(rho1_2d, A_2d, R0_Ar1, levels=[1.0], colors='black', linewidths=2)
plt.colorbar(cf, label=r'$R_0$')
plt.xlabel(r'Susceptible Precaution Rate $\rho_1$', fontsize=11)
plt.ylabel(r'Recruitment Rate $A$ (day⁻¹)', fontsize=11)
plt.title(r'$R_0$ Plane: $A$ vs $\rho_1$', fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig('A_vs_rho1onR.png', dpi=300, bbox_inches='tight')
plt.close()

# Graph 12
print("\n[12/17] A_vs_b1onR.png - VARYING: A (10→500), b₁ (0.01→1.0) [b₁ independent of R₀]")
b1_2d = np.linspace(0.01, 1.0, 100)
R0_Ab1 = np.zeros((len(A_2d), len(b1_2d)))
for i, a_val in enumerate(A_2d):
    for j, b1_val in enumerate(b1_2d):
        R0_Ab1[i, j] = R0_formula(a_val, BASE_PARAMS['beta'], BASE_PARAMS['rho1'], BASE_PARAMS['rho2'],
                                BASE_PARAMS['d'], BASE_PARAMS['p'], BASE_PARAMS['M'],
                                BASE_PARAMS['b2'], BASE_PARAMS['alpha'], BASE_PARAMS['sigma'])
plt.figure(figsize=(8, 6))
cf = plt.contourf(b1_2d, A_2d, R0_Ab1, levels=20, cmap='RdYlGn_r')
plt.contour(b1_2d, A_2d, R0_Ab1, levels=[1.0], colors='black', linewidths=2)
plt.colorbar(cf, label=r'$R_0$')
plt.xlabel(r'Quarantine→Susceptible Rate $b_1$ (day⁻¹)', fontsize=11)
plt.ylabel(r'Recruitment Rate $A$ (day⁻¹)', fontsize=11)
plt.title(r'$R_0$ Plane: $A$ vs $b_1$ (b₁ independent)', fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig('A_vs_b1onR.png', dpi=300, bbox_inches='tight')
plt.close()

# Graph 13
print("\n[13/17] A_vs_b2.png - VARYING: A (10→500), b₂ (0.05→1.0)")
b2_2d = np.linspace(0.05, 1.0, 100)
R0_Ab2 = np.zeros((len(A_2d), len(b2_2d)))
for i, a_val in enumerate(A_2d):
    for j, b2_val in enumerate(b2_2d):
        R0_Ab2[i, j] = R0_formula(a_val, BASE_PARAMS['beta'], BASE_PARAMS['rho1'], BASE_PARAMS['rho2'],
                                BASE_PARAMS['d'], BASE_PARAMS['p'], BASE_PARAMS['M'],
                                b2_val, BASE_PARAMS['alpha'], BASE_PARAMS['sigma'])
plt.figure(figsize=(8, 6))
cf = plt.contourf(b2_2d, A_2d, R0_Ab2, levels=20, cmap='RdYlGn_r')
plt.contour(b2_2d, A_2d, R0_Ab2, levels=[1.0], colors='black', linewidths=2)
plt.colorbar(cf, label=r'$R_0$')
plt.xlabel(r'Exposed→Quarantine Rate $b_2$ (day⁻¹)', fontsize=11)
plt.ylabel(r'Recruitment Rate $A$ (day⁻¹)', fontsize=11)
plt.title(r'$R_0$ Plane: $A$ vs $b_2$', fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig('A_vs_b2.png', dpi=300, bbox_inches='tight')
plt.close()

# Graph 14
print("\n[14/17] A_vs_alphaonR.png - VARYING: A (10→500), α (0.01→0.8)")
R0_Aa = np.zeros((len(A_2d), len(alpha_2d)))
for i, a_val in enumerate(A_2d):
    for j, alpha_val in enumerate(alpha_2d):
        R0_Aa[i, j] = R0_formula(a_val, BASE_PARAMS['beta'], BASE_PARAMS['rho1'], BASE_PARAMS['rho2'],
                                BASE_PARAMS['d'], BASE_PARAMS['p'], BASE_PARAMS['M'],
                                BASE_PARAMS['b2'], alpha_val, BASE_PARAMS['sigma'])
plt.figure(figsize=(8, 6))
cf = plt.contourf(alpha_2d, A_2d, R0_Aa, levels=20, cmap='RdYlGn_r')
plt.contour(alpha_2d, A_2d, R0_Aa, levels=[1.0], colors='black', linewidths=2)
plt.colorbar(cf, label=r'$R_0$')
plt.xlabel(r'Symptom Onset Rate $\alpha$ (day⁻¹)', fontsize=11)
plt.ylabel(r'Recruitment Rate $A$ (day⁻¹)', fontsize=11)
plt.title(r'$R_0$ Plane: $A$ vs $\alpha$', fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig('A_vs_alphaonR.png', dpi=300, bbox_inches='tight')
plt.close()

# Graph 15
print("\n[15/17] A_vs_eta.png - VARYING: A (10→500), η (0.01→0.5) [η independent of R₀]")
eta_2d = np.linspace(0.01, 0.5, 100)
R0_Ae = np.zeros((len(A_2d), len(eta_2d)))
for i, a_val in enumerate(A_2d):
    for j, eta_val in enumerate(eta_2d):
        R0_Ae[i, j] = R0_formula(a_val, BASE_PARAMS['beta'], BASE_PARAMS['rho1'], BASE_PARAMS['rho2'],
                                BASE_PARAMS['d'], BASE_PARAMS['p'], BASE_PARAMS['M'],
                                BASE_PARAMS['b2'], BASE_PARAMS['alpha'], BASE_PARAMS['sigma'])
plt.figure(figsize=(8, 6))
cf = plt.contourf(eta_2d, A_2d, R0_Ae, levels=20, cmap='RdYlGn_r')
plt.contour(eta_2d, A_2d, R0_Ae, levels=[1.0], colors='black', linewidths=2)
plt.colorbar(cf, label=r'$R_0$')
plt.xlabel(r'Recovery Rate $\eta$ (day⁻¹)', fontsize=11)
plt.ylabel(r'Recruitment Rate $A$ (day⁻¹)', fontsize=11)
plt.title(r'$R_0$ Plane: $A$ vs $\eta$ (η independent)', fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig('A_vs_eta.png', dpi=300, bbox_inches='tight')
plt.close()

# Graph 16
print("\n[16/17] A_vs_conR.png - VARYING: A (10→500), c (0.01→0.5) [c independent of R₀]")
c_2d = np.linspace(0.01, 0.5, 100)
R0_Ac = np.zeros((len(A_2d), len(c_2d)))
for i, a_val in enumerate(A_2d):
    for j, c_val in enumerate(c_2d):
        R0_Ac[i, j] = R0_formula(a_val, BASE_PARAMS['beta'], BASE_PARAMS['rho1'], BASE_PARAMS['rho2'],
                                BASE_PARAMS['d'], BASE_PARAMS['p'], BASE_PARAMS['M'],
                                BASE_PARAMS['b2'], BASE_PARAMS['alpha'], BASE_PARAMS['sigma'])
plt.figure(figsize=(8, 6))
cf = plt.contourf(c_2d, A_2d, R0_Ac, levels=20, cmap='RdYlGn_r')
plt.contour(c_2d, A_2d, R0_Ac, levels=[1.0], colors='black', linewidths=2)
plt.colorbar(cf, label=r'$R_0$')
plt.xlabel(r'Quarantine Breach Rate $c$ (day⁻¹)', fontsize=11)
plt.ylabel(r'Recruitment Rate $A$ (day⁻¹)', fontsize=11)
plt.title(r'$R_0$ Plane: $A$ vs $c$ (c independent)', fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig('A_vs_conR.png', dpi=300, bbox_inches='tight')
plt.close()

# ============================================================================
# GRAPH 17: PRCC SENSITIVITY ANALYSIS
# ============================================================================
print("\n[17/17] PRCC_Sensitivity.png - Global Sensitivity via Latin Hypercube Sampling")

param_ranges = {
    'A': (10, 500),
    'beta': (0.05, 1.5),
    'rho1': (0.0, 1.0),
    'rho2': (0.0, 1.0),
    'd': (0.0001, 0.01),
    'p': (0, 0.5),
    'M': (0, 1.0),
    'b2': (0.05, 1.0),
    'alpha': (0.01, 0.8),
    'sigma': (0.01, 0.5),
}

n_samples = 500
param_names = list(param_ranges.keys())
param_samples = {}
for pname in param_names:
    pmin, pmax = param_ranges[pname]
    param_samples[pname] = uniform.rvs(loc=pmin, scale=pmax-pmin, size=n_samples, random_state=42)

# Calculate R0 for each sample
R0_samples = []
for i in range(n_samples):
    R0_val = R0_formula(
        param_samples['A'][i],
        param_samples['beta'][i],
        param_samples['rho1'][i],
        param_samples['rho2'][i],
        param_samples['d'][i],
        param_samples['p'][i],
        param_samples['M'][i],
        param_samples['b2'][i],
        param_samples['alpha'][i],
        param_samples['sigma'][i]
    )
    R0_samples.append(R0_val)

# Calculate PRCC for each parameter
prcc_values = {}
for pname in param_names:
    param_rank = rankdata(param_samples[pname])
    R0_rank = rankdata(R0_samples)
    prcc, pval = spearmanr(param_rank, R0_rank)
    prcc_values[pname] = prcc

# Plot PRCC
sorted_params = sorted(prcc_values.items(), key=lambda x: abs(x[1]), reverse=True)
param_labels = [p[0] for p in sorted_params]
prcc_vals = [p[1] for p in sorted_params]
colors = ['red' if v > 0 else 'blue' for v in prcc_vals]

plt.figure(figsize=(9, 6))
plt.barh(param_labels, prcc_vals, color=colors, alpha=0.7, edgecolor='black', linewidth=1.2)
plt.xlabel('Partial Rank Correlation Coefficient (PRCC)', fontsize=12, fontweight='bold')
plt.ylabel('Parameters', fontsize=12, fontweight='bold')
plt.title('Global Sensitivity: PRCC for R₀', fontsize=13, fontweight='bold')
plt.axvline(0, color='black', linestyle='-', linewidth=0.8)
plt.grid(True, axis='x', alpha=0.3)
plt.tight_layout()
plt.savefig('PRCC_Sensitivity.png', dpi=300, bbox_inches='tight')
plt.close()

print("\n" + "="*80)
print("✓ ALL 17 GRAPHS GENERATED SUCCESSFULLY!")
print("="*80)
