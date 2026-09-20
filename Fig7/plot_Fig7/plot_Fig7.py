import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# Random-matrix reference distribution
# ============================================================
def zd_pdf(k, beta):
    """
    Zakrzewski-Delande curvature distribution

        P(k) = C_beta (1 + k^2)^(-(beta+2)/2),

    where beta = 1 corresponds to the GOE case.
    """
    from math import gamma, sqrt, pi
    C = gamma((beta+2)/2) / (sqrt(pi) * gamma((beta+1)/2))
    return C * (1.0 + k**2) ** (-(beta+2.0)/2.0)

# ============================================================
# Plot formatting
# ============================================================

label_size = 13
tick_size = 11

# ============================================================
# Read processed molecular data
# ============================================================

# Number of electronic states included in the Fig. 7 analysis.
nstates = 200

# Numerical level velocities generated independently for electric
# fields along the x, y, and z directions.
vx = np.loadtxt("v_x_" + str(nstates) + ".txt")
vy = np.loadtxt("v_y_" + str(nstates) + ".txt")
vz = np.loadtxt("v_z_" + str(nstates) + ".txt")

# Numerical level curvatures generated from the same three
# electric-field trajectories.
kappax = np.loadtxt("kappa_x" + str(nstates) + ".txt")
kappay = np.loadtxt("kappa_y" + str(nstates) + ".txt")
kappaz = np.loadtxt("kappa_z" + str(nstates) + ".txt")

# ============================================================
# Normalize level velocities
# ============================================================

# Normalize each Cartesian velocity distribution independently
# so that
#
#     <v_alpha^2> = 1,
#
# with alpha = x, y, z.
Cx = np.nanmean(vx**2)
Cy = np.nanmean(vy**2)
Cz = np.nanmean(vz**2)

v_x = vx / np.sqrt(Cx)
v_y = vy / np.sqrt(Cy)
v_z = vz / np.sqrt(Cz)


# ============================================================
# Combine and normalize curvature data
# ============================================================

# Average the three Cartesian curvature arrays.
#
# This gives one combined molecular curvature dataset representing
# the response to electric-field perturbations along x, y, and z.
k_avg = (kappax + kappay + kappaz) / 3.0

# Normalize the combined curvature distribution by its RMS value:
#
#     K_norm = K / sqrt(<K^2>)
#
# so that the distribution shape can be compared directly with
# the universal random-matrix prediction.
K_rms = np.sqrt(np.nanmean(k_avg**2))
k_norm = k_avg / K_rms


# ============================================================
# Prepare curvature-tail data
# ============================================================

# Use the absolute curvature values for the asymptotic tail analysis.
k_abs = np.abs(k_norm).ravel()

# Remove NaN/Inf entries and zero values before taking logarithms.
k_abs = k_abs[np.isfinite(k_abs)]
k_abs = k_abs[k_abs > 0]

print("Maximum normalized curvature =", k_abs.max())


# Logarithmic bins are used because the curvature distribution
# spans several orders of magnitude.
bins = np.logspace(
    np.log10(k_abs.min()),
    np.log10(k_abs.max()),
    300
)

hist, edges = np.histogram(
    k_abs,
    bins=bins,
    density=True
)

# Use geometric bin centers for logarithmically spaced bins.
centers = np.sqrt(edges[:-1] * edges[1:])

# Remove empty histogram bins before taking logarithms.
mask = hist > 0

lnK = np.log(centers[mask])
lnP = np.log(hist[mask])


# GOE asymptotic curvature tail:
#
#     P(|K|) ~ (2/pi) |K|^-3
#
# or equivalently
#
#     ln P = ln(2/pi) - 3 ln|K|.
lnP_tail = np.log(2.0 / np.pi) - 3.0 * lnK


# ============================================================
# Prepare velocity histogram data
# ============================================================

vals_x = v_x[np.isfinite(v_x)].ravel()
vals_y = v_y[np.isfinite(v_y)].ravel()
vals_z = v_z[np.isfinite(v_z)].ravel()

# Standard normal distribution used as the GOE reference
# for normalized level velocities.
v_grid = np.linspace(-5, 5, 400)

gaussian = np.exp(-v_grid**2 / 2) / np.sqrt(2 * np.pi)


# ============================================================
# Figure 7(a-c): level-velocity distributions
# ============================================================

fig1, axes1 = plt.subplots(
    1,
    3,
    figsize=(10, 3),
    sharey=True
)

ax1, ax2, ax3 = axes1


# ------------------------------------------------------------
# Electric field along x
# ------------------------------------------------------------

ax1.hist(
    vals_x,
    bins=80,
    density=True,
    edgecolor="black"
)

ax1.plot(
    v_grid,
    gaussian,
    color="cyan",
    lw=2,
    label="Gaussian"
)

ax1.set_xlabel(
    r"$v_x$",
    fontsize=label_size
)

ax1.set_ylabel(
    r"$P(v)$",
    fontsize=label_size
)

ax1.set_xlim(-5, 5)

ax1.tick_params(
    axis="both",
    labelsize=tick_size
)

ax1.legend()


# ------------------------------------------------------------
# Electric field along y
# ------------------------------------------------------------

ax2.hist(
    vals_y,
    bins=80,
    density=True,
    edgecolor="black"
)

ax2.plot(
    v_grid,
    gaussian,
    color="cyan",
    lw=2
)

ax2.set_xlabel(
    r"$v_y$",
    fontsize=label_size
)

ax2.set_xlim(-5, 5)

ax2.tick_params(
    axis="both",
    labelsize=tick_size
)

ax2.tick_params(labelleft=False)


# ------------------------------------------------------------
# Electric field along z
# ------------------------------------------------------------

ax3.hist(
    vals_z,
    bins=80,
    density=True,
    edgecolor="black"
)

ax3.plot(
    v_grid,
    gaussian,
    color="cyan",
    lw=2
)

ax3.set_xlabel(
    r"$v_z$",
    fontsize=label_size
)

ax3.set_xlim(-5, 5)

ax3.tick_params(
    axis="both",
    labelsize=tick_size
)

ax3.tick_params(labelleft=False)

fig1.tight_layout()


# ============================================================
# Figure 7(d-f): curvature statistics
# ============================================================

fig2, axes2 = plt.subplots(
    1,
    3,
    figsize=(10, 3)
)

ax4, ax5, ax6 = axes2


# ------------------------------------------------------------
# Figure 7(d): curvature distribution
# ------------------------------------------------------------

vals_c = k_norm[np.isfinite(k_norm)].ravel()

grid = np.linspace(
    np.percentile(vals_c, 0.5),
    np.percentile(vals_c, 99.5),
    800
)

ax4.hist(
    vals_c,
    bins=500,
    density=True,
    alpha=0.4,
    edgecolor="black"
)

# Compare with the GOE Zakrzewski-Delande curvature distribution.
ax4.plot(
    grid,
    zd_pdf(grid, beta=1),
    label=r"$P_{\rm GOE}(K)$"
)

ax4.set_xlabel(
    r"$K$",
    fontsize=label_size
)

ax4.set_ylabel(
    r"$P(K)$",
    fontsize=label_size
)

ax4.set_xlim(-5, 5)

ax4.tick_params(
    axis="both",
    labelsize=tick_size
)

ax4.legend()


# ------------------------------------------------------------
# Figure 7(e): curvature-tail scaling
# ------------------------------------------------------------

ax5.plot(
    lnK,
    lnP,
    "o",
    markersize=3
)

ax5.plot(
    lnK,
    lnP_tail,
    "-",
    linewidth=1,
    color="black",
    label=(
        r"$\ln P = "
        r"\ln\!\left(\frac{2}{\pi}\right)"
        r" - 3\ln|K|$"
    )
)

ax5.set_xlabel(
    r"$\ln|K|$",
    fontsize=label_size
)

ax5.set_ylabel(
    r"$\ln P(|K|)$",
    fontsize=label_size
)

ax5.set_ylim(-10, 10)

ax5.tick_params(
    axis="both",
    labelsize=tick_size
)

ax5.legend()


# ------------------------------------------------------------
# Figure 7(f): GOE-GUE curvature crossover
# ------------------------------------------------------------

# This file is generated independently by plt_Fig7f.py.
# It contains the crossover parameter B in column 1 and the
# corresponding curvature statistic in column 2.
data = np.loadtxt(
    "curvature_vs_B.txt",
    skiprows=1
)

B = data[:, 0]
K2_mean = data[:, 1]

ax6.plot(
    B,
    K2_mean,
    linestyle="-",
    color="black"
)

ax6.set_xlabel(
    r"$B$",
    fontsize=label_size
)

ax6.set_ylabel(
    r"$\langle K^2\rangle$",
    fontsize=label_size
)

# Tick labels are omitted in the published schematic/crossover panel.
ax6.set_xticks([])
ax6.set_yticks([])

fig2.tight_layout()

plt.show()

