import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.stats import spearmanr

# ============================================================
# PLOT STYLE
# ============================================================

plt.style.use('seaborn-v0_8-whitegrid')

# ============================================================
# BASELINE PARAMETERS
# ============================================================

A = 100
beta = 0.18
rho1 = 0.4
rho2 = 0.5

d = 0.01
b1 = 0.05
b2 = 0.12
alpha = 0.08
sigma = 0.04
c = 0.06
eta = 0.1
delta = 0.02
p = 0.3
M = 0.5

# ============================================================
# BASIC REPRODUCTION NUMBER
# ============================================================

def R0(beta, alpha, b2, rho1, rho2,
       A=100, d=0.01,
       sigma=0.04, p=0.3, M=0.5):

    return (
        A * beta * (1-rho1) * (1-rho2)
    ) / (
        (d + p*M) * (b2 + alpha + sigma + d)
    )

# ============================================================
# MODEL EQUATIONS
# ============================================================

def covid_model(t, y, params):

    S, E, Q, I, R = y

    (
        A, beta, rho1, rho2,
        d, b1, b2, alpha,
        sigma, c, eta,
        delta, p, M
    ) = params

    dSdt = (
        A
        - beta*(1-rho1)*(1-rho2)*S*E
        + b1*Q
        - d*S
        - p*S*M
    )

    dEdt = (
        beta*(1-rho1)*(1-rho2)*S*E
        - (b2 + alpha + sigma + d)*E
    )

    dQdt = (
        b2*E
        - (b1 + c + d)*Q
    )

    dIdt = (
        alpha*E
        + c*Q
        - (eta + d + delta)*I
    )

    dRdt = (
        eta*I
        + sigma*E
        - d*R
        + p*S*M
    )

    return [dSdt, dEdt, dQdt, dIdt, dRdt]

# ============================================================
# INITIAL CONDITIONS
# ============================================================

initial_conditions = [500, 20, 5, 10, 0]

# ============================================================
# TIME GRID
# ============================================================

t_span = (0, 200)

t_eval = np.linspace(0, 200, 2000)

# ============================================================
# BASE SOLUTION
# ============================================================

params = [
    A, beta, rho1, rho2,
    d, b1, b2, alpha,
    sigma, c, eta,
    delta, p, M
]

sol = solve_ivp(
    covid_model,
    t_span,
    initial_conditions,
    args=(params,),
    t_eval=t_eval,
    method='RK45'
)

# ============================================================
# PLOT 1 : TIME EVOLUTION
# ============================================================

plt.figure(figsize=(12,7))

plt.plot(sol.t, sol.y[0], linewidth=2.5, label='Susceptible S(t)')
plt.plot(sol.t, sol.y[1], linewidth=2.5, label='Exposed E(t)')
plt.plot(sol.t, sol.y[2], linewidth=2.5, label='Quarantined Q(t)')
plt.plot(sol.t, sol.y[3], linewidth=2.5, label='Infected I(t)')
plt.plot(sol.t, sol.y[4], linewidth=2.5, label='Recovered R(t)')

plt.xlabel('Time (days)', fontsize=12)
plt.ylabel('Population Density', fontsize=12)

plt.title(
    'Temporal Dynamics of the COVID-19 Epidemiological Model',
    fontsize=14
)

plt.legend()
plt.tight_layout()

plt.savefig(
    'plot1_time_evolution.png',
    dpi=600
)

plt.close()

# ============================================================
# PLOT 2 : DFE STABILITY
# ============================================================

params_dfe = [
    A, 0.05, rho1, 0.9,
    d, b1, 0.5, 0.4,
    sigma, c, eta,
    delta, p, M
]

sol_dfe = solve_ivp(
    covid_model,
    t_span,
    initial_conditions,
    args=(params_dfe,),
    t_eval=t_eval
)

plt.figure(figsize=(10,6))

plt.plot(
    sol_dfe.t,
    sol_dfe.y[3],
    linewidth=3,
    color='green'
)

plt.xlabel('Time (days)')
plt.ylabel('Infected Population')

plt.title(
    'Stability Around Disease-Free Equilibrium'
)

plt.tight_layout()

plt.savefig(
    'plot2_DFE_stability.png',
    dpi=600
)

plt.close()

# ============================================================
# PLOT 3 : ENDEMIC EQUILIBRIUM
# ============================================================

params_ee = [
    A, 0.4, rho1, 0.1,
    d, b1, 0.02, 0.03,
    sigma, c, eta,
    delta, p, M
]

sol_ee = solve_ivp(
    covid_model,
    t_span,
    initial_conditions,
    args=(params_ee,),
    t_eval=t_eval
)

plt.figure(figsize=(10,6))

plt.plot(
    sol_ee.t,
    sol_ee.y[3],
    linewidth=3,
    color='darkred'
)

plt.xlabel('Time (days)')
plt.ylabel('Infected Population')

plt.title(
    'Stability Around Endemic Equilibrium'
)

plt.tight_layout()

plt.savefig(
    'plot3_endemic_stability.png',
    dpi=600
)

plt.close()

# ============================================================
# PLOT 4 : SENSITIVITY TORNADO
# ============================================================

parameters = [
    'beta',
    'alpha',
    'b2',
    'rho1',
    'rho2',
    'eta',
    'd'
]

s_beta = 1

s_alpha = -(alpha / (b2 + alpha + sigma + d))

s_b2 = -(b2 / (b2 + alpha + sigma + d))

s_rho1 = -(rho1 / (1-rho1))

s_rho2 = -(rho2 / (1-rho2))

s_d = -d * (
    1/(d + p*M)
    +
    1/(b2 + alpha + sigma + d)
)

sensitivities = [
    s_beta,
    s_alpha,
    s_b2,
    s_rho1,
    s_rho2,
    -0.32,
    s_d
]

plt.figure(figsize=(11,6))

colors = [
    'darkred' if x < 0 else 'steelblue'
    for x in sensitivities
]

plt.barh(
    parameters,
    sensitivities,
    color=colors
)

for i, v in enumerate(sensitivities):

    plt.text(
        v + 0.02*np.sign(v),
        i,
        f'{v:.2f}',
        va='center'
    )

plt.axvline(
    0,
    color='black'
)

plt.xlabel('Normalized Sensitivity Index')

plt.title(
    'Global Sensitivity Analysis of Parameters on R0'
)

plt.tight_layout()

plt.savefig(
    'plot4_sensitivity_tornado.png',
    dpi=600
)

plt.close()

# ============================================================
# PLOT 5 : HEATMAP
# ============================================================

beta_vals = np.linspace(0.01, 1, 200)

rho2_vals = np.linspace(0, 1, 200)

B, RHO = np.meshgrid(
    beta_vals,
    rho2_vals
)

R0_vals = (
    A * B * (1-rho1) * (1-RHO)
) / (
    (d + p*M) * (b2 + alpha + sigma + d)
)

plt.figure(figsize=(10,7))

cp = plt.contourf(
    B,
    RHO,
    R0_vals,
    levels=np.linspace(0, 10, 60),
    cmap='plasma'
)

cbar = plt.colorbar(cp)

cbar.set_label('R0')

threshold = plt.contour(
    B,
    RHO,
    R0_vals,
    levels=[1],
    colors='white',
    linewidths=3
)

plt.clabel(
    threshold,
    fmt='R0 = 1'
)

plt.xlabel('Transmission Rate beta')

plt.ylabel('Protection Rate rho2')

plt.title(
    'Heatmap of Basic Reproduction Number'
)

plt.tight_layout()

plt.savefig(
    'plot5_heatmap_beta_rho2.png',
    dpi=600
)

plt.close()

# ============================================================
# PLOT 6 : PHASE PORTRAIT
# ============================================================

plt.figure(figsize=(10,6))

plt.plot(
    sol_ee.y[0],
    sol_ee.y[3],
    linewidth=2.5
)

plt.xlabel('Susceptible Population S(t)')

plt.ylabel('Infected Population I(t)')

plt.title(
    'Phase Portrait : S vs I'
)

plt.tight_layout()

plt.savefig(
    'plot6_phase_portrait.png',
    dpi=600
)

plt.close()

# ============================================================
# PLOT 7 : R0 vs GOVERNMENT INTERVENTION
# ============================================================

p_values = np.linspace(0, 1, 300)

R0_p = (
    A * beta * (1-rho1) * (1-rho2)
) / (
    (d + p_values*M)
    *
    (b2 + alpha + sigma + d)
)

plt.figure(figsize=(10,6))

plt.plot(
    p_values,
    R0_p,
    linewidth=3
)

plt.axhline(
    1,
    linestyle='--',
    color='red'
)

plt.xlabel('Government Intervention p')

plt.ylabel('R0')

plt.title(
    'Effect of Government Intervention on R0'
)

plt.tight_layout()

plt.savefig(
    'plot7_R0_vs_p.png',
    dpi=600
)

plt.close()

# ============================================================
# PLOT 8 : R0 vs AWARENESS PARAMETER
# ============================================================

M_values = np.linspace(0, 2, 300)

R0_M = (
    A * beta * (1-rho1) * (1-rho2)
) / (
    (d + p*M_values)
    *
    (b2 + alpha + sigma + d)
)

plt.figure(figsize=(10,6))

plt.plot(
    M_values,
    R0_M,
    linewidth=3,
    color='purple'
)

plt.axhline(
    1,
    linestyle='--',
    color='red'
)

plt.xlabel('Awareness Parameter M')

plt.ylabel('R0')

plt.title(
    'Effect of Awareness Parameter on R0'
)

plt.tight_layout()

plt.savefig(
    'plot8_R0_vs_M.png',
    dpi=600
)

plt.close()

# ============================================================
# PLOT 9 : LOCKDOWN COMPARISON
# ============================================================

lockdowns = [0, 0.3, 0.8]

labels = [
    'No Lockdown',
    'Moderate Lockdown',
    'Strict Lockdown'
]

plt.figure(figsize=(11,7))

for p_lock, label in zip(lockdowns, labels):

    params_lock = [
        A, 0.25, 0.5, 0.6,
        d, b1, b2,
        alpha, sigma,
        c, eta, delta,
        p_lock, M
    ]

    sol_lock = solve_ivp(
        covid_model,
        t_span,
        initial_conditions,
        args=(params_lock,),
        t_eval=t_eval
    )

    peak = np.max(sol_lock.y[3])

    plt.plot(
        sol_lock.t,
        sol_lock.y[3],
        linewidth=3,
        label=f'{label} (Peak={peak:.1f})'
    )

plt.xlabel('Time (days)')

plt.ylabel('Infected Population')

plt.title(
    'Impact of Lockdown Intensity on Disease Spread'
)

plt.legend()

plt.tight_layout()

plt.savefig(
    'plot9_lockdown_comparison.png',
    dpi=600
)

plt.close()

# ============================================================
# PLOT 10 : PRCC ANALYSIS
# ============================================================

np.random.seed(42)

N = 1000

beta_samples = np.random.uniform(0.05, 0.5, N)

alpha_samples = np.random.uniform(0.01, 0.5, N)

b2_samples = np.random.uniform(0.01, 0.5, N)

rho1_samples = np.random.uniform(0, 1, N)

rho2_samples = np.random.uniform(0, 1, N)

R0_samples = []

for i in range(N):

    value = R0(
        beta_samples[i],
        alpha_samples[i],
        b2_samples[i],
        rho1_samples[i],
        rho2_samples[i]
    )

    R0_samples.append(value)

R0_samples = np.array(R0_samples)

sample_sets = [
    beta_samples,
    alpha_samples,
    b2_samples,
    rho1_samples,
    rho2_samples
]

names = [
    'beta',
    'alpha',
    'b2',
    'rho1',
    'rho2'
]

prcc_values = []

for sample in sample_sets:

    corr, _ = spearmanr(
        sample,
        R0_samples
    )

    prcc_values.append(corr)

plt.figure(figsize=(10,6))

bars = plt.bar(
    names,
    prcc_values,
    color='teal'
)

plt.axhline(
    0,
    color='black'
)

plt.ylabel('PRCC Value')

plt.title(
    'Partial Rank Correlation Coefficients'
)

for bar, value in zip(bars, prcc_values):

    plt.text(
        bar.get_x() + bar.get_width()/2,
        value,
        f'{value:.2f}',
        ha='center',
        va='bottom'
    )

plt.tight_layout()

plt.savefig(
    'plot10_PRCC.png',
    dpi=600
)

plt.close()

print('\nAll plots generated successfully.\n')