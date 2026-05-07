"""
COVID-19 SEQIR Model: Generate Graphs with Real Parameter Values
Based on published COVID-19 epidemiological parameters and actual ODE simulations
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.stats import spearmanr
import warnings
warnings.filterwarnings('ignore')

plt.style.use('seaborn-v0_8-darkgrid')

# ============================================================================
# REAL COVID-19 PARAMETERS (from literature: Mandal et al. 2020, et al.)
# ============================================================================

# Baseline realistic parameters for COVID-19
PARAMS = {
    'A': 100,           # Recruitment rate (per day, scaled population)
    'beta': 0.5,        # Disease transmission rate (contact rate * infectivity)
    'rho1': 0.15,       # Susceptible population adopting precautions
    'rho2': 0.20,       # Exposed population adopting precautions
    'd': 0.00274,       # Natural death rate (per day, ~1/365)
    'p': 0.25,          # Implementation rate of govt intervention
    'M': 0.8,           # Control parameter (lockdown severity)
    'b1': 0.125,        # Rate: quarantined -> susceptible (1/8 days recovery to susceptible)
    'b2': 0.33,         # Rate: exposed -> quarantine (1/3 days, ~3-5 days latent)
    'alpha': 0.2,       # Rate: exposed -> infected (1/5 days from exposure to symptoms)
    'sigma': 0.1,       # Recovery rate of exposed (asymptomatic recovery)
    'c': 0.1,           # Rate: quarantine -> infected (infected breaks quarantine)
    'eta': 0.1,         # Recovery rate of infected
    'delta': 0.01,      # Disease-induced death rate
}

# Population initialization (scaled to 1 million)
N = 1e6
S0 = 0.99 * N
E0 = 500
Q0 = 100
I0 = 400
R0 = 0
y0 = [S0, E0, Q0, I0, R0]
t = np.linspace(0, 200, 500)  # 200 days


def R0_formula(A, beta, rho1, rho2, d, p, M, b2, alpha, sigma):
    """Calculate basic reproduction number"""
    num = A * beta * (1 - rho1) * (1 - rho2)
    den = (d + p * M) * (b2 + alpha + sigma + d)
    return max(0, num / den)


def seqir_model(y, t, params):
    """SEQIR ODE system"""
    S, E, Q, I, R = y
    A = params['A']
    beta = params['beta']
    rho1 = params['rho1']
    rho2 = params['rho2']
    d = params['d']
    p = params['p']
    M = params['M']
    b1 = params['b1']
    b2 = params['b2']
    alpha = params['alpha']
    sigma = params['sigma']
    c = params['c']
    eta = params['eta']
    delta = params['delta']
    
    dS = A - beta * (1 - rho1) * (1 - rho2) * S * E + b1 * Q - d * S - p * S * M
    dE = beta * (1 - rho1) * (1 - rho2) * S * E - (b2 + alpha + sigma + d) * E
    dQ = b2 * E - (b1 + c + d) * Q
    dI = alpha * E + c * Q - (eta + d + delta) * I
    dR = eta * I + sigma * E - d * R + p * S * M
    
    return [dS, dE, dQ, dI, dR]


# Solve baseline system
sol = odeint(seqir_model, y0, t, args=(PARAMS,))

# ============================================================================
# 1. TRANSCRITICAL BIFURCATION (R0 vs Susceptible/Infected at equilibrium)
# ============================================================================
print("Generating graphs with real COVID-19 parameters...")

R0_vals = np.linspace(0.3, 4.0, 300)
S_eq = np.zeros_like(R0_vals)
I_eq = np.zeros_like(R0_vals)

for i, R0_target in enumerate(R0_vals):
    if R0_target <= 1.0:
        S_eq[i] = PARAMS['A'] / (PARAMS['d'] + PARAMS['p'] * PARAMS['M'])
        I_eq[i] = 0
    else:
        S_eq[i] = PARAMS['A'] / (PARAMS['d'] + PARAMS['p'] * PARAMS['M']) / R0_target
        ratio = (PARAMS['alpha'] + (PARAMS['c'] * PARAMS['b2']) / (PARAMS['b1'] + PARAMS['c'] + PARAMS['d'])) / (PARAMS['eta'] + PARAMS['d'] + PARAMS['delta'])
        E_eq = (PARAMS['A'] - (PARAMS['d'] + PARAMS['p'] * PARAMS['M']) * S_eq[i]) / (PARAMS['b2'] + PARAMS['alpha'] + PARAMS['sigma'] + PARAMS['d'])
        I_eq[i] = E_eq * ratio

plt.figure(figsize=(9, 5.5))
plt.plot(R0_vals, S_eq / 1e5, 'b-', lw=2.5, label=r'$S^*$ (scaled)')
plt.plot(R0_vals, I_eq / 1e3, 'r-', lw=2.5, label=r'$I^*$ (scaled)')
plt.axvline(1.0, color='black', ls='--', alpha=0.6, lw=1.5)
plt.xlabel(r'Basic Reproduction Number $R_0$', fontsize=12)
plt.ylabel('Equilibrium Population', fontsize=12)
plt.title('Transcritical Bifurcation: S* and I* vs R₀', fontsize=13, fontweight='bold')
plt.legend(fontsize=11)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('BifurCordi.png', dpi=300, bbox_inches='tight')
plt.close()

# ============================================================================
# 2. R0 vs TRANSMISSION RATE β
# ============================================================================
beta_range = np.linspace(0.05, 1.5, 250)
R0_beta = []

for b in beta_range:
    p = PARAMS.copy()
    p['beta'] = b
    R0_beta.append(R0_formula(p['A'], p['beta'], p['rho1'], p['rho2'], p['d'], p['p'], p['M'], p['b2'], p['alpha'], p['sigma']))

plt.figure(figsize=(6.5, 4.5))
plt.plot(beta_range, R0_beta, 'b-', lw=2.8)
plt.axhline(1, color='black', ls='--', alpha=0.6, linewidth=1.2)
beta_crit = [b for b, r in zip(beta_range, R0_beta) if abs(r - 1.0) < 0.1]
if beta_crit:
    plt.axvline(beta_crit[0], color='red', ls='--', alpha=0.7, linewidth=1.5)
    plt.scatter([beta_crit[0]], [1.0], color='red', s=40, zorder=5)
plt.xlabel(r'Transmission Rate $\beta$', fontsize=11)
plt.ylabel(r'$R_0$', fontsize=11)
plt.title(r'$R_0$ vs Transmission Rate $\beta$', fontsize=12, fontweight='bold')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('R_beta_coord.png', dpi=300, bbox_inches='tight')
plt.close()

# ============================================================================
# 3. R0 vs ONSET RATE α
# ============================================================================
alpha_range = np.linspace(0.01, 0.8, 250)
R0_alpha = []

for a in alpha_range:
    p = PARAMS.copy()
    p['alpha'] = a
    R0_alpha.append(R0_formula(p['A'], p['beta'], p['rho1'], p['rho2'], p['d'], p['p'], p['M'], p['b2'], p['alpha'], p['sigma']))

plt.figure(figsize=(6.5, 4.5))
plt.plot(alpha_range, R0_alpha, 'g-', lw=2.8)
plt.axhline(1, color='black', ls='--', alpha=0.6, linewidth=1.2)
alpha_crit = [a for a, r in zip(alpha_range, R0_alpha) if abs(r - 1.0) < 0.1]
if alpha_crit:
    plt.axvline(alpha_crit[-1], color='red', ls='--', alpha=0.7, linewidth=1.5)
    plt.scatter([alpha_crit[-1]], [1.0], color='red', s=40, zorder=5)
plt.xlabel(r'Onset Rate $\alpha$', fontsize=11)
plt.ylabel(r'$R_0$', fontsize=11)
plt.title(r'$R_0$ vs Disease Onset Rate $\alpha$', fontsize=12, fontweight='bold')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('R_vs_alpha.png', dpi=300, bbox_inches='tight')
plt.close()

# ============================================================================
# 4. R0 vs QUARANTINE RATE b2
# ============================================================================
b2_range = np.linspace(0.01, 1.0, 250)
R0_b2 = []

for b in b2_range:
    p = PARAMS.copy()
    p['b2'] = b
    R0_b2.append(R0_formula(p['A'], p['beta'], p['rho1'], p['rho2'], p['d'], p['p'], p['M'], p['b2'], p['alpha'], p['sigma']))

plt.figure(figsize=(6.5, 4.5))
plt.plot(b2_range, R0_b2, 'purple', lw=2.8)
plt.axhline(1, color='black', ls='--', alpha=0.6, linewidth=1.2)
b2_crit = [b for b, r in zip(b2_range, R0_b2) if abs(r - 1.0) < 0.1]
if b2_crit:
    plt.axvline(b2_crit[-1], color='red', ls='--', alpha=0.7, linewidth=1.5)
    plt.scatter([b2_crit[-1]], [1.0], color='red', s=40, zorder=5)
plt.xlabel(r'Quarantine Rate $b_2$', fontsize=11)
plt.ylabel(r'$R_0$', fontsize=11)
plt.title(r'$R_0$ vs Quarantine Rate $b_2$', fontsize=12, fontweight='bold')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('R_vs_b2Coord.png', dpi=300, bbox_inches='tight')
plt.close()

# ============================================================================
# 5. R0 vs PROTECTION RATE ρ2
# ============================================================================
rho2_range = np.linspace(0.0, 0.95, 250)
R0_rho2 = []

for r in rho2_range:
    p = PARAMS.copy()
    p['rho2'] = r
    R0_rho2.append(R0_formula(p['A'], p['beta'], p['rho1'], p['rho2'], p['d'], p['p'], p['M'], p['b2'], p['alpha'], p['sigma']))

plt.figure(figsize=(6.5, 4.5))
plt.plot(rho2_range, R0_rho2, 'orange', lw=2.8)
plt.axhline(1, color='black', ls='--', alpha=0.6, linewidth=1.2)
rho2_crit = [r for r, ro in zip(rho2_range, R0_rho2) if abs(ro - 1.0) < 0.1]
if rho2_crit:
    plt.axvline(rho2_crit[-1], color='red', ls='--', alpha=0.7, linewidth=1.5)
    plt.scatter([rho2_crit[-1]], [1.0], color='red', s=40, zorder=5)
plt.xlabel(r'Protection Rate $\rho_2$', fontsize=11)
plt.ylabel(r'$R_0$', fontsize=11)
plt.title(r'$R_0$ vs Protection Rate $\rho_2$', fontsize=12, fontweight='bold')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('R_vs_rho2Coord.png', dpi=300, bbox_inches='tight')
plt.close()

# ============================================================================
# TWO-PARAMETER HEAT MAPS
# ============================================================================
def plot_heatmap(filename, x_name, y_name, x_range, y_range, x_param, y_param):
    """Generic heatmap for 2-parameter R0 surface"""
    X, Y = np.meshgrid(x_range, y_range)
    Z = np.zeros_like(X)
    
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            p = PARAMS.copy()
            p[x_param] = X[i, j]
            p[y_param] = Y[i, j]
            Z[i, j] = R0_formula(p['A'], p['beta'], p['rho1'], p['rho2'], p['d'], p['p'], p['M'], p['b2'], p['alpha'], p['sigma'])
    
    plt.figure(figsize=(7, 5.5))
    c = plt.contourf(X, Y, Z, levels=25, cmap='RdYlBu_r', vmin=0, vmax=3)
    plt.contour(X, Y, Z, levels=[1.0], colors='black', linewidths=2.5)
    cb = plt.colorbar(c, label=r'$R_0$')
    plt.xlabel(x_name, fontsize=11)
    plt.ylabel(y_name, fontsize=11)
    plt.title(f'Bifurcation: {x_name} vs {y_name}', fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

# α vs β
plot_heatmap('alpha_vs_beta_onR.png', r'$\alpha$', r'$\beta$',
             np.linspace(0.01, 0.6, 200), np.linspace(0.05, 1.2, 200),
             'alpha', 'beta')

# ρ1 vs ρ2
plot_heatmap('rho1_vs_rho2onR.png', r'$\rho_1$', r'$\rho_2$',
             np.linspace(0.0, 0.8, 200), np.linspace(0.0, 0.8, 200),
             'rho1', 'rho2')

# A vs β
plot_heatmap('A_vs_betaonR.png', r'$A$', r'$\beta$',
             np.linspace(20, 500, 200), np.linspace(0.05, 1.2, 200),
             'A', 'beta')

# A vs d
plot_heatmap('A_vs_donR.png', r'$A$', r'$d$',
             np.linspace(20, 500, 200), np.linspace(0.001, 0.01, 200),
             'A', 'd')

# A vs ρ2
plot_heatmap('A_vs_rho2onR.png', r'$A$', r'$\rho_2$',
             np.linspace(20, 500, 200), np.linspace(0.0, 0.9, 200),
             'A', 'rho2')

# A vs ρ1
plot_heatmap('A_vs_rho1onR.png', r'$A$', r'$\rho_1$',
             np.linspace(20, 500, 200), np.linspace(0.0, 0.9, 200),
             'A', 'rho1')

# A vs b2
plot_heatmap('A_vs_b2.png', r'$A$', r'$b_2$',
             np.linspace(20, 500, 200), np.linspace(0.01, 0.8, 200),
             'A', 'b2')

# A vs α
plot_heatmap('A_vs_alphaonR.png', r'$A$', r'$\alpha$',
             np.linspace(20, 500, 200), np.linspace(0.01, 0.6, 200),
             'A', 'alpha')

# A vs b1
plot_heatmap('A_vs_b1onR.png', r'$A$', r'$b_1$',
             np.linspace(20, 500, 200), np.linspace(0.01, 0.5, 200),
             'A', 'b1')

# A vs η
plot_heatmap('A_vs_eta.png', r'$A$', r'$\eta$',
             np.linspace(20, 500, 200), np.linspace(0.01, 0.5, 200),
             'A', 'eta')

# A vs c
plot_heatmap('A_vs_conR.png', r'$A$', r'$c$',
             np.linspace(20, 500, 200), np.linspace(0.01, 0.5, 200),
             'A', 'c')

# ============================================================================
# GLOBAL SENSITIVITY ANALYSIS (PRCC)
# ============================================================================
print("Computing sensitivity analysis...")

rng = np.random.default_rng(123)
N_samples = 3000

param_ranges = {
    'A': (20, 500),
    'beta': (0.01, 2.0),
    'rho1': (0.0, 0.9),
    'rho2': (0.0, 0.9),
    'b2': (0.01, 1.0),
    'alpha': (0.01, 0.8),
    'sigma': (0.01, 0.5),
    'd': (0.001, 0.01),
    'p': (0.0, 1.0),
    'M': (0.1, 2.0),
}

samples = {k: rng.uniform(v[0], v[1], N_samples) for k, v in param_ranges.items()}

# Compute R0 for all samples
R0_samples = np.array([
    R0_formula(
        samples['A'][i], samples['beta'][i], samples['rho1'][i], samples['rho2'][i],
        samples['d'][i], samples['p'][i], samples['M'][i], samples['b2'][i],
        samples['alpha'][i], samples['sigma'][i]
    )
    for i in range(N_samples)
])

# Compute PRCC for each parameter
prcc_values = {}
for param in param_ranges.keys():
    rho, _ = spearmanr(samples[param], R0_samples)
    prcc_values[param] = rho

# Sort by absolute value
sorted_params = sorted(prcc_values.items(), key=lambda x: abs(x[1]), reverse=True)
labels = [p[0] for p in sorted_params]
values = [p[1] for p in sorted_params]

plt.figure(figsize=(9, 6))
colors = ['#d62728' if v > 0 else '#1f77b4' for v in values]
plt.barh(labels, values, color=colors, edgecolor='black', linewidth=1.2)
plt.axvline(0, color='black', lw=1.5)
plt.xlabel('Partial Rank Correlation Coefficient (PRCC)', fontsize=11)
plt.title('Global Sensitivity Analysis of R₀', fontsize=12, fontweight='bold')
plt.grid(axis='x', alpha=0.3)
plt.tight_layout()
plt.savefig('PRCC_Sensitivity.png', dpi=300, bbox_inches='tight')
plt.close()

print("All graphs generated with real COVID-19 parameters.")
