#python3
import re
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# Gap-ratio statistics
# ============================================================
def gap_ratios(delta):
    """
    Compute adjacent-gap ratios

        r_n = min(delta_n, delta_{n+1}) / max(delta_n, delta_{n+1})

    from a sequence of nearest-neighbor level spacings.

    Unlike the spacing distribution P(s), the gap-ratio statistic does
    not require spectral unfolding.
    """

    delta = np.asarray(delta)

    with np.errstate(divide='ignore', invalid='ignore'):
        r = np.minimum(delta[1:], delta[:-1]) / \
            np.maximum(delta[1:], delta[:-1])

    # Ignore undefined ratios arising from exactly degenerate levels (0/0).
    r_mean = np.nanmean(r)

    # Fraction of ratios that are NaN
    nan_fraction = np.mean(np.isnan(r))

    print("Mean gap ratio <r> =", r_mean)
    print("NaN fraction =", nan_fraction)

    return r, r_mean, nan_fraction


# ============================================================
# Unfolding diagnostics
# ============================================================

def plot_unfold_fit(E, N, N_smooth, title="Unfolding fit N(E)"):
    """Compare the exact spectral staircase with its smooth polynomial fit."""
    plt.figure()
    plt.plot(E, N, '.', label='data: (E_i, i)')
    plt.plot(E, N_smooth, '-', label='smooth fit N_smooth(E)')
    plt.xlabel("E")
    plt.ylabel("N(E)")
    plt.title(title)
    plt.legend()
    plt.show()

def plot_spacings(xi_sorted, title="Unfolded spacings"):
    """Plot the nearest-neighbor spacings of the unfolded spectrum."""
    s = np.diff(xi_sorted)
    plt.figure()
    plt.hist(s, bins=50)
    plt.xlabel("s = xi_{i+1} - xi_i")
    plt.ylabel("count")
    plt.title(title)
    plt.show()

    print("min(s) =", s.min(), "mean(s) =", s.mean(), "std(s) =", s.std())

# ============================================================
# Random-matrix reference distributions
# ============================================================
def wigner_dyson(s, B):
    return 2 * B * s * np.exp(-B * s**2)

def poisson_DOS(s):
    return np.exp(-s)

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

    # Diagnostic plots used to inspect the quality of the unfolding. 
    plot_spacings(xi_sorted)
    plot_unfold_fit(E, N, xi_sorted)
    return xi

# ============================================================
# Density of states
# ============================================================
def DOS_plots(data,cut_states0,cut_states1,bin_width,filename):
    """
    Calculate and plot the density of states of the original spectrum.

    The vertical dashed lines indicate the spectral window retained
    for the level-spacing analysis. The histogram data are written
    to <filename>_DOS.txt and are used to construct Supplemental
    Figure S1.
    """     

    data = np.array(data)
    bins = np.arange(data.min(), data.max() + bin_width, bin_width)
    hist_y, hist_x = np.histogram(data, bins=bins, density=True)
    x_centers = 0.5 * (hist_x[:-1] + hist_x[1:])
    
    hist_area = np.sum(hist_y * bin_width)
    print("Histogram area =", hist_area)

    fig, axes = plt.subplots(1, 1, figsize=(8,5))
    axes.bar(x_centers, hist_y, width=bin_width, align='center', edgecolor='k')
   
    # Mark the lower and upper boundaries of the spectrum used
    # for the level-spacing analysis.
    axes.axvline(x=data[cut_states0], color='r', linestyle='--')
    axes.axvline(x=data[-cut_states1], color='r', linestyle='--')

    axes.set_xlabel(r"Orbital Energy $\epsilon$ (Hartree)")
    axes.set_ylabel(r"Probability density $P(\epsilon)$")

    axes.legend()
    fig.text(0.5, 0.98, "Density of States", 
     ha='center', va='center', fontsize=10)
 
    # Save the DOS histogram used by the final plotting script.   
    np.savetxt(
     filename+"_DOS.txt",
     np.column_stack((x_centers, hist_y)),
     delimiter=",",
     header=str(data[cut_states0])+"_"+str(data[-cut_states1]),
     comments=""
    )

# ============================================================
# Nearest-neighbor level-spacing statistics
# ============================================================
def wignor_dyson_plot(data,data_unfold,bin_width1,bin_width2,filename):
    """
    Calculate nearest-neighbor level spacings before and after
    spectral unfolding.

    The unfolded spacing distribution is compared with the
    Poisson and random-matrix reference distributions.

    The histogram of unfolded spacings is written to
    <filename>.txt. These processed data are used by the final
    plotting script for Figure 1.
    """

    # Nearest-neighbor spacings before and after unfolding.
    spacings = np.diff(data)
    spacings_unfold = np.diff(data_unfold)

    # Normalize to unit mean spacing.
    s = spacings / np.mean(spacings)
    s_unfold = spacings_unfold / np.mean(spacings_unfold)
 
    print("s mean {} ".format(s.mean()))
    print("s_unfold mean {} ".format(s_unfold.mean()))
   
    smin = min(spacings)
    smin_unfold = min(spacings_unfold)
    print("smin {} smin_unfold {} ".format(smin,smin_unfold))

    # Gap ratios are reported for both the raw spectra.
    # The raw-spectrum value provides an unfolding-independent
    # diagnostic of the spectral statistics.
    gap_ratios(spacings)

    bins_1 = np.arange(min(spacings), max(spacings) + bin_width1, bin_width1)
    bins_2 = np.arange(min(spacings_unfold), max(spacings_unfold) + bin_width2, bin_width2)
    hist_y1, hist_x1 = np.histogram(spacings, bins='auto', density=True)
    hist_y2, hist_x2 = np.histogram(spacings_unfold, bins=bins_2, density=True)
    x_centers1 = 0.5 * (hist_x1[:-1] + hist_x1[1:])
    x_centers2 = 0.5 * (hist_x2[:-1] + hist_x2[1:])

    # Reference spacing distributions.
    s_plot = np.linspace(0, 10, 400)
    P_GOE = wigner_dyson(s_plot, np.pi/4.0)
    P_POS = poisson_DOS(s_plot)

    # Diagnostic comparison before and after unfolding.
    fig, axes = plt.subplots(2, 1, figsize=(8,5))
    widths = np.diff(hist_x1)
    axes[0].bar(x_centers1, hist_y1,width=widths, align='center', edgecolor='k')
    axes[0].set_ylabel("Probability density")
    axes[0].legend()
    axes[0].set_title("Before unfolding")
    axes[1].bar(x_centers2, hist_y2, width=bin_width2, align='center', edgecolor='k')
    axes[1].plot(s_plot,P_GOE,color='r',label="Wigner–Dyson (GOE)")
    axes[1].plot(s_plot,P_POS,color='orange' ,label="Poisson")
    axes[1].set_xlabel("Unfolded level spacing s = ΔE / ⟨ΔE⟩")
    axes[1].set_ylabel("Probability density p(s)")
    axes[1].set_xlim(-1,5)
    axes[1].legend()
    axes[1].set_title("After unfolding")
    fig.text(0.5, 0.98, "Neighboring Eigenvalue Difference", 
     ha='center', va='center', fontsize=10)
  
    # Save the processed spacing histogram used to construct Fig. 1. 
    np.savetxt(
     filename+".txt",
     np.column_stack((x_centers2, hist_y2)),
     delimiter=",",
     header="x_center,hist_y",
     comments=""
    )
 


# ============================================================
# Input files and analysis parameters
# ============================================================

# Select the molecular system to analyze.
# <filename>.dat contains the CIS excited state energies
# extracted from the Q-Chem calculation.
# <filename>.param contains the parameters used in the analysis.
#filename = "18_7N_dz"
filename = "18_v6_5_dz"

with open(filename+".dat","r") as f:
    text=f.read()
numbers = re.findall(r"[-+]?\d*\.\d+|\d+", text)
data = np.array(numbers, dtype=float)


# Read analysis parameters.
#
# Entries in <filename>.param:
#   param[0] : polynomial degree used for spectral unfolding
#   param[1] : fraction of lowest orbital-energy levels excluded 
#   param[2] : fraction of highest orbital-energy levels excluded 
param = np.loadtxt(filename+".param")
order = int(param[0])
cut0 = int(param[1])
cut1 = float(param[2])

# If plotting only occupied orbitals
# read in number of occupied orbitals 
Nocc = 0
plt_nocc = 0
if (len(param)==4):
   plt_nocc = 1
   Nocc = int(param[3]) 


N = len(data)
data = sorted(data)

# Determine the orbital-energy window used for the single-particle
# level-spacing analysis. cut0 and cut1 specify the fractions removed
# from the low- and high-energy ends of the orbital spectrum.
cut_states0 = int(cut0*N)
cut_states1 = int(cut1*N)

# Keep the complete orbital spectrum for the DOS plot.
raw_data = data

# Retain only the selected spectral window for unfolding and
# nearest-neighbor level-spacing analysis.
data= data[cut_states0:-cut_states1]

print("cut {} states from the lower end and {} from the upper end".format(cut_states0,cut_states1))
print("Length of data before cut {} after cut {} ".format(len(raw_data),len(data)))

# Unfold the selected orbital-energy spectrum using the polynomial
# unfolding procedure.
data_unfolded = unfolding(data,order)

# ------------------------------------------------------------
# Density of states
# ------------------------------------------------------------

# Generate the DOS of the complete orbital spectrum.
# The selected spectral window is indicated in the DOS plot.
# The corresponding histogram data are saved for Supplemental Fig. S1.
bin_width = 0.002
DOS_plots(raw_data,cut_states0,cut_states1,bin_width,filename)

# ------------------------------------------------------------
# Nearest-neighbor level-spacing statistics
# ------------------------------------------------------------

# Histogram bin widths before and after unfolding.
bin_width_1 = 0.01
bin_width_2 = 0.08

# Calculate the nearest-neighbor spacing statistics and save the
# processed unfolded-spacing histogram used for Fig. 1.
wignor_dyson_plot(data,data_unfolded,bin_width_1,bin_width_2,filename)

# Calculate the nearest-neighbor spacing statistics for the occupied orbitals and save the
# processed unfolded-spacing histogram used for Fig. 8.
bin_width_1 = 0.0005
bin_width_2 = 0.05
if (plt_nocc == 1):
    data = data[0:Nocc]
    data_unfolded = data_unfolded[0:Nocc]
    wignor_dyson_plot(data,data_unfolded,bin_width_1,bin_width_2,filename+"_nocc")

plt.show()

