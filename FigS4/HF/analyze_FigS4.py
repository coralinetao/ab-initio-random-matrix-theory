import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# Gap-ratio statistic
# ============================================================
def gap_ratios(delta):
    """
    Compute the mean adjacent-gap ratio

        r_n = min(delta_n, delta_{n+1})
              / max(delta_n, delta_{n+1})

    from a sequence of nearest-neighbor level spacings.

    The gap-ratio statistic is useful because it does not require
    spectral unfolding.
    """

    delta = np.asarray(delta)

    with np.errstate(divide="ignore", invalid="ignore"):
        r = (
            np.minimum(delta[1:], delta[:-1])
            / np.maximum(delta[1:], delta[:-1])
        )

    # Ignore undefined ratios arising from exactly degenerate levels.
    r_mean = np.nanmean(r)

    print("Mean gap ratio <r> =", r_mean)

    return r_mean

# ============================================================
# Random-matrix reference distributions
# ============================================================

def wigner_dyson_goe(s):
    """
    GOE Wigner-Dyson nearest-neighbor spacing distribution.
    """
    return (np.pi / 2.0) * s * np.exp(-np.pi * s**2 / 4.0)

def poisson_spacing(s):
    """
    Poisson nearest-neighbor spacing distribution.
    """
    return np.exp(-s)


# ============================================================
# Polynomial spectral unfolding
# ============================================================

def unfolding(data, degree):
    """
    Unfold a spectrum using a polynomial fit to the integrated
    density of states.

    For ordered energies E_i, define the spectral staircase

        N(E_i) = i.

    A smooth polynomial approximation to N(E) is fitted and used to
    define the unfolded energies

        xi_i = N_smooth(E_i).

    The energy axis is centered and scaled before fitting to improve
    numerical conditioning.
    """

    E = np.sort(data)
    N = np.arange(1, len(E) + 1)

    # Center and scale the energy variable.
    E0 = E.mean()
    Escale = E.std() if E.std() > 0 else 1.0
    x = (E - E0) / Escale

    # Polynomial fit to the smooth spectral staircase.
    coefficients = np.polyfit(x, N, degree)

    # Unfolded energies.
    xi = np.polyval(coefficients, x)

    # Diagnostic: unfolded spacings should remain positive.
    spacings = np.diff(xi)

    print("Minimum unfolded spacing =", np.min(spacings))
    print("Mean unfolded spacing =", np.mean(spacings))
    print("Number of non-positive spacings =", np.sum(spacings <= 0))

    return xi


# ============================================================
# Nearest-neighbor spacing analysis
# ============================================================

def spacing_statistics(data, data_unfolded, bin_width, filename):
    """
    Calculate nearest-neighbor level-spacing statistics before and
    after unfolding.

    The unfolded spacing histogram is saved to <filename>.txt.
    """

    spacings = np.diff(data)
    spacings_unfolded = np.diff(data_unfolded)

    # Normalize unfolded spacings to unit mean.
    s_unfolded = (
        spacings_unfolded
        / np.mean(spacings_unfolded)
    )

    # Gap ratio from the original spectrum.
    # This statistic is independent of the unfolding procedure.
    r_mean = gap_ratios(spacings)

    # Histogram of unfolded spacings.
    bins = np.arange(
        s_unfolded.min(),
        s_unfolded.max() + bin_width,
        bin_width
    )

    hist_y, hist_x = np.histogram(
        s_unfolded,
        bins=bins,
        density=True
    )

    x_centers = 0.5 * (
        hist_x[:-1] + hist_x[1:]
    )

    # Reference distributions.
    s_plot = np.linspace(0, 5, 400)

    P_GOE = wigner_dyson_goe(s_plot)
    P_POS = poisson_spacing(s_plot)

    # Diagnostic plot.
    plt.figure()

    plt.bar(
        x_centers,
        hist_y,
        width=bin_width,
        align="center",
        edgecolor="k"
    )

    plt.plot(
        s_plot,
        P_GOE,
        label="GOE"
    )

    plt.plot(
        s_plot,
        P_POS,
        label="Poisson"
    )

    plt.xlabel("Unfolded level spacing s")
    plt.ylabel("Probability density")
    plt.xlim(0, 5)
    plt.legend()

    # Save processed spacing data.
    np.savetxt(
        filename + ".txt",
        np.column_stack(
            (x_centers, hist_y)
        ),
        delimiter=",",
        header="x_center,hist_y",
        comments=""
    )

    return r_mean


# ============================================================
# Analysis parameters
# ============================================================

# All ten orbital spectra use the same analysis parameters.
#
# The parameter file contains:
#
#   param[0] : polynomial degree used for unfolding
#   param[1] : fraction of lowest levels excluded
#   param[2] : fraction of highest levels excluded
#   param[3] : optional number of occupied levels retained
#
param = np.loadtxt(
    "benzene_perturbed_acqz.param"
)

degree = int(param[0])
cut0 = float(param[1])
cut1 = float(param[2])

# Optional restriction to the first Nocc levels.
use_nocc = len(param) == 4

if use_nocc:
    Nocc = int(param[3])


# ============================================================
# Analyze orbital-resolved spectra
# ============================================================

ravg = []

for i in range(10):

    filename = "orbital" + str(i + 1)

    # Read the orbital-resolved spectrum.
    data = np.loadtxt(
        filename + ".dat"
    )

    # If requested, retain only the specified number of levels.
    if use_nocc:
        data = data[:Nocc]

    data = np.sort(data)

    N = len(data)

    # Remove the low- and high-energy spectral edges before
    # performing the level-statistics analysis.
    cut_states0 = int(cut0 * N)
    cut_states1 = int(cut1 * N)

    data_selected = data[
        cut_states0:-cut_states1
    ]

    print(
        "\n{}: {} levels before cut, {} after cut".format(
            filename,
            len(data),
            len(data_selected)
        )
    )

    # Unfold the selected spectrum.
    data_unfolded = unfolding(
        data_selected,
        degree
    )

    # Calculate spacing statistics and save the processed histogram.
    r_mean = spacing_statistics(
        data_selected,
        data_unfolded,
        bin_width=0.05,
        filename=filename
    )

    # Save the unfolded spectrum for reproducibility.
    np.savetxt(
        filename + "_unfolded.dat",
        data_unfolded
    )

    ravg.append(r_mean)


# ============================================================
# Figure S4: average gap ratio vs orbital index
# ============================================================

orbital_index = np.arange(
    1,
    len(ravg) + 1
)

plt.figure()

plt.plot(
    orbital_index,
    ravg,
    "o-"
)

plt.xlabel("Orbital")
plt.ylabel("Average gap ratio")

plt.tight_layout()
plt.show()
