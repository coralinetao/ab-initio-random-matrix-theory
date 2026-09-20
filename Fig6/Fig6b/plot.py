import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

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

    s = np.diff(xi)
    print("Minimum unfolded spacing =", np.min(s))

    xi_sorted = xi[seq]

    return xi

# ============================================================
# Analysis parameters
# ============================================================

# Polynomial degree used for unfolding.
degree = 13

# State window displayed in Fig. 6(b).
#
# Python slicing excludes the upper endpoint, so this selects
# states 1230 through 1269.
state_min = 480
state_max = 521

# Parameter frames included in the calculation.
frame_min = 1
frame_max = 250

# Increment of the external parameter between consecutive frames,
# in atomic units.
dlam = 2.0 

# ============================================================
# Read parameter-dependent electronic energies
# ============================================================

# E.dat is organized as
#
#     rows    -> electronic states
#     columns -> time steps
#
E_states_x_frames = np.loadtxt("E.dat")

nstates_tot = E_states_x_frames.shape[0]

print("Shape of E.dat =", E_states_x_frames.shape)

# ============================================================
# Select the spectral window used for unfolding
# ============================================================

# The unfolding is performed using the central 60% of the complete
# electronic spectrum.
#
# This window is much larger than the 40 states displayed in Fig. 6(b)
# so that the smooth variation of the density of states can be
# determined from a broader spectral region.
nunfold_1 = int(0.25 * nstates_tot)
nunfold_2 = int(0.70 * nstates_tot)

print(
    "Unfolding using states {} through {}".format(
        nunfold_1,
        nunfold_2 - 1
    )
)

# Raw energies in the spectral region used for unfolding.
E_raw = E_states_x_frames[
    nunfold_1:nunfold_2,
    frame_min:frame_max
]

nframes = E_raw.shape[1]


# Allocate the unfolded spectrum.
E_unfolded = np.full_like(E_raw, np.nan)

# ============================================================
# Unfold the spectrum independently at each parameter value
# ============================================================
for step in range(nframes):
    Ei_sorted = E_raw[:,step]

    E_unfolded[:,step] = unfolding(Ei_sorted,degree)

# ============================================================
# Extract the states displayed in Fig. 6(b)
# ============================================================

# Convert the absolute state indices into indices relative to the
# larger unfolding window.
slice_start = state_min - nunfold_1
slice_stop = state_max - nunfold_1

# Unfolded energies of the states displayed in the figure.
E_unfolded_slice = E_unfolded[
    slice_start:slice_stop,
    :
]

# Corresponding raw electronic energies.
E_slice = E_raw[
    slice_start:slice_stop,
    :
]

# ============================================================
# Construct the external-parameter coordinate
# ============================================================

# Original parameter values:
#
#     lambda = 0, dlam, 2*dlam, ...
#
lam = np.arange(nframes, dtype=float) * dlam


# ============================================================
# Calculate level velocities
# ============================================================

# Time coordinate corresponding to the columns of E.dat.
time = np.arange(nframes, dtype=float) * dlam

# Level velocity of the unfolded spectrum:
#
#     V_i(t) = d epsilon_i / dt
#
# The derivative is evaluated using a Savitzky-Golay filter.
V = savgol_filter(
    E_unfolded_slice,
    window_length=7,
    polyorder=3,
    deriv=1,
    delta=dlam,
    axis=1
)


# ============================================================
# Remove boundary points
# ============================================================

# Exclude points near the boundaries where the numerical derivative
# is less reliable.
trim = 5

time_trim = time[trim:-trim]
E_trim = E_unfolded_slice[:, trim:-trim]
V_trim = V[:, trim:-trim]


# ============================================================
# Construct the rescaled parameter x
# ============================================================

# Mean-square level velocity:
#
#     C(0) = <V_i(t)^2>
#
# where the average is over the selected states and time points.
C0 = np.nanmean(V_trim**2)

# Characteristic level-velocity scale.
s = np.sqrt(C0)

# Rescaled parameter used on the x-axis of Fig. 6(b):
#
#     x = sqrt(C(0)) (t - t_0)
#
x = s * (time_trim - time_trim[0])

print("C(0) =", C0)
print("sqrt(C(0)) =", s)


# ============================================================
# Plot Fig. 6(b)
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
for i in range(E_unfolded_slice.shape[0]):   # loop over states
      axes[0].plot(lam, E_slice[i,:], lw=1,color='black')
      axes[1].plot(x, E_unfolded_slice[i,trim:-trim], lw=1,color='black')

plt.xlabel("Rescaled Parameter x", fontsize=13)
plt.ylabel("Unfolded energy levels ε", fontsize=13)
for ax in axes.flat:
    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", labelsize=14)
plt.show()

