import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# Polynomial spectral unfolding
# ============================================================
def unfolding(data,degree):
    """
    Unfold an energy spectrum using a polynomial fit to the
    integrated density of states (spectral staircase).

    For ordered energies E_i, define

        N(E_i) = i.

    A polynomial is fitted to N(E), giving the smooth staircase
    N_smooth(E). The unfolded energies are

        xi_i = N_smooth(E_i).

    The energy axis is centered and scaled before fitting to improve
    numerical conditioning of the polynomial fit.
    """

    E = np.sort(data)
    N = np.arange(1, len(E)+1)
    seq = np.argsort(E)

    # Center and scale the energy variable for numerical stability.
    E0 = E.mean()
    Escale = E.std() if E.std() > 0 else 1.0
    x = (E - E0) / Escale

    # Fit the spectral staircase N(E).
    c = np.polyfit(x, N, degree)
    xi = np.polyval(c, x)

    # Smooth DOS:
    #
    # rho_smooth(E) = dN_smooth(E)/dE
    #
    # Because the polynomial is fitted as a function of
    # x = (E-E0)/Escale, the factor 1/Escale follows from
    # the chain rule.
    dc = np.polyder(c)
    rho = np.polyval(dc, x) / Escale

    # Local smooth mean level spacing:
    #
    # Delta(E) = 1 / rho_smooth(E)
    Delta = 1.0 / rho             # local mean spacing in ENERGY units
    Delta0 = np.median(Delta)
    Delta_mean = np.mean(Delta)
    Delta_std  = np.std(Delta)
    rel_var = Delta_std / Delta_mean

    print("Single energy scale median {} mean {} sd {} sd/mean {}".format(Delta0,Delta_mean,Delta_std,rel_var))
    s = np.diff(xi)

    xi_sorted = xi[seq]

    return xi


# ============================================================
# Input parameters
# ============================================================

# State window used for the Fig. 7 and Fig. S5 analysis.
#
# Python slicing uses i0:i1, so state i1 is excluded.
i0 = 425
i1 = 625

# Polynomial degree used for spectral unfolding.
degree = 8

# ============================================================
# Read electronic energies
# ============================================================

# E.dat:
#
#     rows    -> electronic states
#     columns -> successive electric-field configurations
#
E_raw = np.loadtxt("E.dat")

nframes = E_raw.shape[1]

print("E.dat shape =", E_raw.shape)
print("Selected states =", i0, "to", i1 - 1)
print("Number of selected states =", i1 - i0)

# ============================================================
# Construct the electric-field path coordinate
# ============================================================

# efield_history.dat contains the electric-field trajectory.
#
# Column 0 gives the step index.
# Columns 1-3 give the Cartesian electric-field components.
efield = np.loadtxt("efield_history.dat")

field = efield[:, 1:4]

# Distance between successive points along the electric-field path:
#
#     d lambda =
#       | E(t+1) - E(t) |
#
# The cumulative path length is used as the parameter lambda.
dfield = np.linalg.norm(
    field[1:] - field[:-1],
    axis=1
)

lam = np.zeros(nframes)

for step in range(1, nframes):
    lam[step] = np.sum(dfield[:step])

# ============================================================
# Unfold the selected spectrum at every field configuration
# ============================================================

# Select the same state window at every frame.
E_window = E_raw[i0:i1, :]

# Allocate array for unfolded energies.
E_unfolded = np.full_like(
    E_window,
    np.nan
)

for step in range(nframes):

    E_unfolded[:, step] = unfolding(
        E_window[:, step],
        degree
    )


# ============================================================
# Numerical level velocities
# ============================================================

# Level velocity with respect to the cumulative electric-field
# path coordinate:
#
#     V_i(lambda) = d epsilon_i / d lambda
#
V_ref = np.gradient(
    E_unfolded,
    lam,
    axis=1,
    edge_order=2
)


# ============================================================
# Numerical level curvatures
# ============================================================

# Level curvature:
#
#     K_i(lambda) = d^2 epsilon_i / d lambda^2
#
K_ref = np.gradient(
    V_ref,
    lam,
    axis=1,
    edge_order=2
)


# ============================================================
# Save processed data
# ============================================================

nstates = i1 - i0

np.savetxt(
    "v_x_" + str(nstates) + ".txt",
    V_ref
)

np.savetxt(
    "kappa_x" + str(nstates) + ".txt",
    K_ref
)

print(
    "Saved v_x_{}.txt and kappa_x{}.txt".format(
        nstates,
        nstates
    )
)


