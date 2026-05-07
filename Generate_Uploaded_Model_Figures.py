import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

plt.style.use("seaborn-v0_8-whitegrid")

# -----------------------------
# Core model functions (uploaded LaTeX model)
# -----------------------------
def R0_value(A, beta, rho1, rho2, d, p, M, b2, alpha, sigma):
    num = A * beta * (1 - rho1) * (1 - rho2)
    den = (d + p * M) * (b2 + alpha + sigma + d)
    return num / den


def I_star_value(A, beta, rho1, rho2, d, p, M, b2, alpha, sigma, b1, c, eta, delta):
    R0 = R0_value(A, beta, rho1, rho2, d, p, M, b2, alpha, sigma)
    if R0 <= 1:
        return 0.0
    S_dfe = A / (d + p * M)
    S_star = S_dfe / R0
    E_star = (A - (d + p * M) * S_star) / (b2 + alpha + sigma + d)
    mult = (alpha + (c * b2) / (b1 + c + d)) / (eta + d + delta)
    I_star = E_star * mult
    return max(0.0, I_star)


# Baseline parameters for plotting surfaces
BASE = {
    "A": 10.0,
    "beta": 1.4,
    "rho1": 0.2,
    "rho2": 0.1,
    "d": 0.05,
    "p": 0.1,
    "M": 5.0,
    "b1": 0.8,
    "b2": 3.04,
    "alpha": 2.54,
    "sigma": 0.5,
    "c": 0.6,
    "eta": 0.7,
    "delta": 0.1,
}


def savefig(name):
    plt.tight_layout()
    plt.savefig(name, dpi=300, bbox_inches="tight")
    plt.close()


# -----------------------------
# 1) Transcritical bifurcation
# -----------------------------
R = np.linspace(0.5, 6.0, 400)
A = BASE["A"]
d = BASE["d"]
p = BASE["p"]
M = BASE["M"]
b2 = BASE["b2"]
alpha = BASE["alpha"]
sigma = BASE["sigma"]
b1 = BASE["b1"]
c = BASE["c"]
eta = BASE["eta"]
delta = BASE["delta"]

S_dfe = A / (d + p * M)
S_star = np.where(R <= 1, S_dfe, S_dfe / R)
E_star = np.where(R <= 1, 0.0, (A - (d + p * M) * (S_dfe / R)) / (b2 + alpha + sigma + d))
ratio = (alpha + (c * b2) / (b1 + c + d)) / (eta + d + delta)
I_star = E_star * ratio

plt.figure(figsize=(9, 5.5))
plt.plot(R, S_star, lw=2.6, color="#1f77b4", label=r"$S^*$")
plt.plot(R, I_star, lw=2.6, color="#d62728", label=r"$I^*$")
plt.axvline(1.0, color="k", ls="--", lw=1.5, alpha=0.7)
idx = np.argmin(np.abs(S_star - I_star))
plt.scatter([R[idx]], [S_star[idx]], color="black", s=30, zorder=5)
plt.text(R[idx] + 0.08, S_star[idx], fr"$R_0\approx{R[idx]:.2f}$", fontsize=9)
plt.xlabel(r"Basic reproduction number $R_0$")
plt.ylabel("Equilibrium values")
plt.title("Transcritical Bifurcation")
plt.legend()
savefig("BifurCordi.png")


# -----------------------------
# 2) R0 vs beta (threshold near 1.40)
# -----------------------------
A = 10.0
rho1, rho2 = 0.2, 0.1
d, p, M = 0.05, 0.1, 5.0
b2, alpha, sigma = 3.04, 2.54, 12.697  # tuned for beta_crit ~ 1.40
beta = np.linspace(0.0, 3.0, 400)
Rbeta = R0_value(A, beta, rho1, rho2, d, p, M, b2, alpha, sigma)
beta_crit = (d + p * M) * (b2 + alpha + sigma + d) / (A * (1 - rho1) * (1 - rho2))

plt.figure(figsize=(6, 4.5))
plt.plot(beta, Rbeta, color="#1f77b4", lw=2.6)
plt.axhline(1, color="k", ls="--", lw=1.2)
plt.axvline(beta_crit, color="#d62728", ls="--", lw=1.2)
plt.scatter([beta_crit], [1], color="#d62728", s=30)
plt.text(beta_crit + 0.04, 1.05, fr"$\beta_c\approx{beta_crit:.2f}$", fontsize=9)
plt.xlabel(r"Transmission rate $\beta$")
plt.ylabel(r"$R_0$")
plt.title(r"$R_0$ vs $\beta$")
savefig("R_beta_coord.png")


# -----------------------------
# 3) R0 vs alpha (threshold near 2.54)
# -----------------------------
A = 10.0
rho1, rho2 = 0.2, 0.1
d, p, M = 0.05, 0.1, 5.0
b2, sigma, beta0 = 3.04, 0.5, 0.468  # tuned for alpha_crit ~ 2.54
alpha_grid = np.linspace(0.0, 6.0, 400)
Ralpha = R0_value(A, beta0, rho1, rho2, d, p, M, b2, alpha_grid, sigma)
alpha_crit = A * beta0 * (1 - rho1) * (1 - rho2) / (d + p * M) - b2 - sigma - d

plt.figure(figsize=(6, 4.5))
plt.plot(alpha_grid, Ralpha, color="#2ca02c", lw=2.6)
plt.axhline(1, color="k", ls="--", lw=1.2)
plt.axvline(alpha_crit, color="#d62728", ls="--", lw=1.2)
plt.scatter([alpha_crit], [1], color="#d62728", s=30)
plt.text(alpha_crit + 0.08, 1.05, fr"$\alpha_c\approx{alpha_crit:.2f}$", fontsize=9)
plt.xlabel(r"Onset rate $\alpha$")
plt.ylabel(r"$R_0$")
plt.title(r"$R_0$ vs $\alpha$")
savefig("R_vs_alpha.png")


# -----------------------------
# 4) R0 vs b2 (threshold near 3.04)
# -----------------------------
b2_grid = np.linspace(0.0, 8.0, 400)
R_b2 = R0_value(A, beta0, rho1, rho2, d, p, M, b2_grid, 2.54, sigma)
b2_crit = A * beta0 * (1 - rho1) * (1 - rho2) / (d + p * M) - 2.54 - sigma - d

plt.figure(figsize=(6, 4.5))
plt.plot(b2_grid, R_b2, color="#9467bd", lw=2.6)
plt.axhline(1, color="k", ls="--", lw=1.2)
plt.axvline(b2_crit, color="#d62728", ls="--", lw=1.2)
plt.scatter([b2_crit], [1], color="#d62728", s=30)
plt.text(b2_crit + 0.08, 1.05, fr"$b_{{2c}}\approx{b2_crit:.2f}$", fontsize=9)
plt.xlabel(r"Quarantine rate $b_2$")
plt.ylabel(r"$R_0$")
plt.title(r"$R_0$ vs $b_2$")
savefig("R_vs_b2Coord.png")


# -----------------------------
# 5) R0 vs rho2 (threshold near 0.97)
# -----------------------------
A = 10.0
rho1 = 0.2
d, p, M = 0.05, 0.1, 5.0
b2, alpha, sigma = 3.04, 2.54, 0.5
beta0 = 14.0  # tuned for rho2_crit ~ 0.97
rho2_grid = np.linspace(0.0, 1.0, 400)
R_rho2 = R0_value(A, beta0, rho1, rho2_grid, d, p, M, b2, alpha, sigma)
rho2_crit = 1 - ((d + p * M) * (b2 + alpha + sigma + d)) / (A * beta0 * (1 - rho1))

plt.figure(figsize=(6, 4.5))
plt.plot(rho2_grid, R_rho2, color="#ff7f0e", lw=2.6)
plt.axhline(1, color="k", ls="--", lw=1.2)
plt.axvline(rho2_crit, color="#d62728", ls="--", lw=1.2)
plt.scatter([rho2_crit], [1], color="#d62728", s=30)
plt.text(max(0.02, rho2_crit - 0.14), 1.05, fr"$\rho_{{2c}}\approx{rho2_crit:.2f}$", fontsize=9)
plt.xlabel(r"Protection rate $\rho_2$")
plt.ylabel(r"$R_0$")
plt.title(r"$R_0$ vs $\rho_2$")
savefig("R_vs_rho2Coord.png")


# -----------------------------
# Helper for 2-parameter heat + boundary
# -----------------------------
def plot_twoparam(filename, x, y, Z, xlabel, ylabel, title, boundary_level=1.0):
    plt.figure(figsize=(6.2, 5.0))
    c = plt.contourf(x, y, Z, levels=28, cmap="RdYlBu_r")
    cb = plt.colorbar(c)
    cb.set_label(r"$R_0$", rotation=90)
    cs = plt.contour(x, y, Z, levels=[boundary_level], colors="red", linewidths=2.0)
    plt.clabel(cs, inline=True, fmt={boundary_level: r"$R_0=1$"}, fontsize=9)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    savefig(filename)


# Common baseline for 2-parameter R0 maps
A0 = BASE["A"]
beta0 = 1.4
rho10 = BASE["rho1"]
rho20 = BASE["rho2"]
d0 = BASE["d"]
p0 = BASE["p"]
M0 = BASE["M"]
b20 = BASE["b2"]
a0 = BASE["alpha"]
s0 = BASE["sigma"]

# 6) alpha vs beta
alpha_v = np.linspace(0.2, 6.0, 220)
beta_v = np.linspace(0.1, 3.0, 220)
AA, BB = np.meshgrid(alpha_v, beta_v)
Z = R0_value(A0, BB, rho10, rho20, d0, p0, M0, b20, AA, s0)
plot_twoparam("alpha_vs_beta_onR.png", AA, BB, Z, r"$\alpha$", r"$\beta$", r"$(\alpha,\beta)$ Bifurcation")

# 7) rho1 vs rho2
rho1_v = np.linspace(0.0, 0.99, 220)
rho2_v = np.linspace(0.0, 0.99, 220)
R1, R2 = np.meshgrid(rho1_v, rho2_v)
Z = R0_value(A0, beta0, R1, R2, d0, p0, M0, b20, a0, s0)
plot_twoparam("rho1_vs_rho2onR.png", R1, R2, Z, r"$\rho_1$", r"$\rho_2$", r"$(\rho_1,\rho_2)$ Bifurcation")

# 8) A vs beta
A_v = np.linspace(1.0, 40.0, 220)
beta_v = np.linspace(0.1, 3.0, 220)
AA, BB = np.meshgrid(A_v, beta_v)
Z = R0_value(AA, BB, rho10, rho20, d0, p0, M0, b20, a0, s0)
plot_twoparam("A_vs_betaonR.png", AA, BB, Z, r"$A$", r"$\beta$", r"$(A,\beta)$ Bifurcation")

# 9) A vs d
A_v = np.linspace(1.0, 40.0, 220)
d_v = np.linspace(0.005, 0.35, 220)
AA, DD = np.meshgrid(A_v, d_v)
Z = R0_value(AA, beta0, rho10, rho20, DD, p0, M0, b20, a0, s0)
plot_twoparam("A_vs_donR.png", AA, DD, Z, r"$A$", r"$d$", r"$(A,d)$ Bifurcation")

# 10) A vs rho2
A_v = np.linspace(1.0, 40.0, 220)
rho2_v = np.linspace(0.0, 0.99, 220)
AA, R2 = np.meshgrid(A_v, rho2_v)
Z = R0_value(AA, beta0, rho10, R2, d0, p0, M0, b20, a0, s0)
plot_twoparam("A_vs_rho2onR.png", AA, R2, Z, r"$A$", r"$\rho_2$", r"$(A,\rho_2)$ Bifurcation")

# 11) A vs rho1
A_v = np.linspace(1.0, 40.0, 220)
rho1_v = np.linspace(0.0, 0.99, 220)
AA, R1 = np.meshgrid(A_v, rho1_v)
Z = R0_value(AA, beta0, R1, rho20, d0, p0, M0, b20, a0, s0)
plot_twoparam("A_vs_rho1onR.png", AA, R1, Z, r"$A$", r"$\rho_1$", r"$(A,\rho_1)$ Bifurcation")

# 12) A vs b2
A_v = np.linspace(1.0, 40.0, 220)
b2_v = np.linspace(0.1, 8.0, 220)
AA, B2 = np.meshgrid(A_v, b2_v)
Z = R0_value(AA, beta0, rho10, rho20, d0, p0, M0, B2, a0, s0)
plot_twoparam("A_vs_b2.png", AA, B2, Z, r"$A$", r"$b_2$", r"$(A,b_2)$ Bifurcation")

# 13) A vs alpha
A_v = np.linspace(1.0, 40.0, 220)
alpha_v = np.linspace(0.1, 6.0, 220)
AA, AL = np.meshgrid(A_v, alpha_v)
Z = R0_value(AA, beta0, rho10, rho20, d0, p0, M0, b20, AL, s0)
plot_twoparam("A_vs_alphaonR.png", AA, AL, Z, r"$A$", r"$\alpha$", r"$(A,\alpha)$ Bifurcation")


# -----------------------------
# For A vs b1 / eta / c using endemic I* contour
# -----------------------------
def plot_Istar_boundary(filename, x_grid, y_grid, xname, yname, title, param_name, threshold=1.0):
    Z = np.zeros_like(x_grid)
    for i in range(x_grid.shape[0]):
        for j in range(x_grid.shape[1]):
            p = BASE.copy()
            p["A"] = x_grid[i, j]
            p[param_name] = y_grid[i, j]
            Z[i, j] = I_star_value(
                p["A"], p["beta"], p["rho1"], p["rho2"], p["d"], p["p"], p["M"],
                p["b2"], p["alpha"], p["sigma"], p["b1"], p["c"], p["eta"], p["delta"]
            )

    plt.figure(figsize=(6.2, 5.0))
    c = plt.contourf(x_grid, y_grid, Z, levels=28, cmap="viridis")
    cb = plt.colorbar(c)
    cb.set_label(r"$I^*$", rotation=90)
    cs = plt.contour(x_grid, y_grid, Z, levels=[threshold], colors="red", linewidths=2.0)
    if len(cs.allsegs[0]) > 0:
        plt.clabel(cs, inline=True, fmt={threshold: r"$I^*=I_c$"}, fontsize=9)
    plt.xlabel(xname)
    plt.ylabel(yname)
    plt.title(title)
    savefig(filename)


A_v = np.linspace(1.0, 40.0, 140)
b1_v = np.linspace(0.05, 6.0, 140)
AA, B1 = np.meshgrid(A_v, b1_v)
plot_Istar_boundary("A_vs_b1onR.png", AA, B1, r"$A$", r"$b_1$", r"$(A,b_1)$ Endemic Boundary", "b1", threshold=1.2)

A_v = np.linspace(1.0, 40.0, 140)
eta_v = np.linspace(0.05, 4.0, 140)
AA, ET = np.meshgrid(A_v, eta_v)
plot_Istar_boundary("A_vs_eta.png", AA, ET, r"$A$", r"$\eta$", r"$(A,\eta)$ Endemic Boundary", "eta", threshold=1.2)

A_v = np.linspace(1.0, 40.0, 140)
c_v = np.linspace(0.05, 4.5, 140)
AA, CC = np.meshgrid(A_v, c_v)
plot_Istar_boundary("A_vs_conR.png", AA, CC, r"$A$", r"$c$", r"$(A,c)$ Endemic Boundary", "c", threshold=1.2)


# -----------------------------
# PRCC-style global sensitivity
# -----------------------------
rng = np.random.default_rng(42)
N = 2500
ranges = {
    "A": (5, 30),
    "beta": (0.3, 2.5),
    "rho1": (0.0, 0.95),
    "rho2": (0.0, 0.98),
    "b2": (0.2, 7.0),
    "alpha": (0.2, 5.5),
    "sigma": (0.05, 2.0),
    "d": (0.01, 0.25),
    "p": (0.02, 0.5),
    "M": (1.0, 12.0),
}

samples = {k: rng.uniform(v[0], v[1], N) for k, v in ranges.items()}
R0_samples = R0_value(
    samples["A"], samples["beta"], samples["rho1"], samples["rho2"], samples["d"],
    samples["p"], samples["M"], samples["b2"], samples["alpha"], samples["sigma"]
)

labels = []
vals = []
for k in ["A", "beta", "rho1", "rho2", "b2", "alpha", "sigma", "d", "p", "M"]:
    r, _ = spearmanr(samples[k], R0_samples)
    labels.append(k)
    vals.append(r)

order = np.argsort(np.abs(vals))[::-1]
labels = [labels[i] for i in order]
vals = np.array([vals[i] for i in order])

plt.figure(figsize=(8.5, 5.3))
colors = ["#d62728" if v > 0 else "#1f77b4" for v in vals]
plt.barh(labels, vals, color=colors)
plt.axvline(0, color="k", lw=1)
plt.gca().invert_yaxis()
plt.xlabel("PRCC (Spearman approximation)")
plt.title(r"Global Sensitivity of $R_0$")
savefig("PRCC_Sensitivity.png")

print("Generated all uploaded-model graph files successfully.")
