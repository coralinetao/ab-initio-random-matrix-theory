#python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter


# ============================================================
# Input unfolded spectra
# ============================================================

# Each .dat file contains the unfolded many-particle energy levels
# for one molecular system. These are the same unfolded CIS spectra
# used for the many-particle level-statistics analysis.
filenames = ["alanine_adz_4000_unfolded","benzene1_adz_B0_unfolded","c6h6_adz_4000_unfolded","ch3_oxi_adz_3000_unfolded","ph_nh2_ch3_adz_unfolded"]

# Labels corresponding to the spectra listed above.
labels = ["Alanine", "Benzene (C1)", "Benzene (D6h)", "Methyloxirane", "Phenylethylamine"]

# ============================================================
# Spectral-form-factor parameters
# ============================================================
# Dimensionless time:
#
#     tau = t / t_H
#
# where t_H is the Heisenberg time.
#
# Because the spectra have been unfolded to unit mean level spacing,
# Delta = 1 and t_H = 2*pi (hbar = 1). The factor 2*pi therefore
# appears explicitly in the Fourier phase below.

t = np.linspace(0, 2.0, 5000)

# Spectral-form-factor results are averaged over overlapping
# windows of the unfolded spectrum to reduce fluctuations.
window_size = 300
step = 50

count = 0
# ============================================================
# Calculate the spectral form factor
# ============================================================
for filename in filenames: 

     # Read and sort the unfolded energy levels.
     data = np.loadtxt(filename+".dat")
     data = np.sort(data)  # always numpy array, always sorted 
     N = len(data)

     # Accumulator for averaging K(tau) over spectral windows.
     K_accum = np.zeros_like(t, dtype=float)
     n_windows = 0

     # Divide the unfolded spectrum into overlapping windows.
     # Each window contains 300 consecutive levels and successive
     # windows are shifted by 50 levels.
     for start in range(0, N - window_size + 1, step):
         w = data[start:start + window_size]
         print("start {} end {}".format(start,start + window_size))

         M = len(w)

        # Spectral Fourier amplitude:
        #
        #     Z(tau) = sum_n exp(-2*pi*i*tau*epsilon_n)
        #
        # where epsilon_n are unfolded energy levels.
         Z = np.exp(-2j * np.pi * t[:, None] * w[None, :]).sum(axis=1)

        # Normalized spectral form factor:
        #
        #     K(tau) = |Z(tau)|^2 / M^2
        #
        # With this normalization K(0) = 1.
         K = (np.abs(Z) ** 2) / (M * M)

         K_accum += K
         n_windows += 1
     print("n_windows {}".format(n_windows))

     # Average the spectral form factor over all windows.
     K_avg = K_accum / n_windows
     plt.plot(t, K_avg, label=labels[count],alpha=1)

     # For the present normalization, the long-time plateau
     # associated with a window containing M levels is 1/M.
     plt.axhline(1.0/M, linestyle='--', linewidth=1)
     count += 1      
   
# ============================================================
# Plot Figure 4
# ============================================================
plt.xlim(-0.04, 2.0) 
plt.xlabel(r"$\tau(t/t_H)$", fontsize=16)  
plt.ylabel(r"$K(\tau)$", fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.yscale("log")
plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.legend(fontsize=14) 
#plt.savefig("spectral_form_factor.png", dpi=300, bbox_inches="tight")  # saves figure
plt.show()
