#python3
import re
import numpy as np
import itertools
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline, PchipInterpolator
from scipy.stats import kstest
from scipy.integrate import quad
from scipy.special import gamma
from itertools import combinations
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def wigner_dyson(s, B):
    return 2 * B * s * np.exp(-B * s**2)

def wigner_dyson_GUE(s):
    return 32.0/np.pi/np.pi* s**2 * np.exp(-4.0/np.pi * s**2)

def wigner_dyson_GSE(s):
    return (2**18)/(3**6 * np.pi**3) * s**4 * np.exp(-64.0/(9*np.pi) * s**2)

def wigner_DOS(s):
    return 0.5 / np.pi*np.sqrt(4.0-s**2) 

def poisson_DOS(s):
    return np.exp(-s)

labels = ["(a)", "(b)", "(c)",
          "(d)", "(e)", "(f)",
          "(g)", "(h)", "(i)"]

files = [["ch3_oxi_acqz.txt","alanine_acqz.txt"],["ph_nh2_ch3_acqz.txt",
"benzene_perturbed_acqz.txt"]]
#files = [["ch3_oxi_adz_1500.txt","alanine_adz_4000.txt"],
#["ph_nh2_ch3_adz.txt","benzene1_adz_B0.txt"]
#]

#file_inset = "c6h6_adz_4000.txt"
file_inset = "benzene_eq_acqz.txt"
bin_width2 = 0.05
s_plot = np.linspace(0, 10, 400)
P_GOE = wigner_dyson(s_plot, np.pi/4.0)
P_POS = poisson_DOS(s_plot)


fig, axes = plt.subplots(2, 2, figsize=(10, 5), sharex=True, sharey=True)

for i in range(2):
    for j in range(2):
        x_centers2, hist_y2 = np.loadtxt(
            files[i][j],
            delimiter=",",
            skiprows=1,
            unpack=True
        )

        ax = axes[i, j]
        ax.bar(x_centers2, hist_y2, width=bin_width2, align='center', edgecolor='k', linewidth=0.1,color='black')
        ax.plot(s_plot,P_GOE,color='r',label="Wigner–Dyson (GOE)")
        ax.plot(s_plot,P_POS,color='orange' ,label="Wigner–Dyson (GUE)")
        ax.set_xlim(-0.2, 6)
        ax.set_ylim(0.0,1.2)


# Load inset data
x_inset, y_inset = np.loadtxt(
    file_inset,
    delimiter=",",
    skiprows=1,
    unpack=True
)

# Select bottom-right subplot
ax_main = axes[1, 1]

# Create inset inside it
ax_inset = inset_axes(ax_main, width="50%", height="60%", loc="upper right")

# Plot inset histogram
ax_inset.bar(x_inset, y_inset,
             width=bin_width2,
             align='center',
             edgecolor='k',
             linewidth=0.1, color = 'black')

# Plot reference curves
ax_inset.plot(s_plot, P_GOE,color='red')
ax_inset.plot(s_plot, P_POS,color='orange')

# Optional zoom (typical for spacing distributions)
ax_inset.set_xlim(-0.2, 4)

# Make ticks smaller
ax_inset.tick_params(labelsize=8)
ax_inset.legend().remove()

ax_inset.text(
    0.05, 0.95, "(e)",
    transform=ax_inset.transAxes,
    fontsize=12,
    verticalalignment='top',
    horizontalalignment='left'
)

fig.supxlabel(r"Unfolded level spacing $s$", fontsize=18)
fig.supylabel(r"Probability density $P(s)$", fontsize=18)



#handles, labels = axes[0,0].get_legend_handles_labels()
#   
#fig.legend(
#    handles, labels,
#    loc="upper center",
#    ncol=2,          # put GOE and Poisson side by side
#    frameon=False
#)  

for ax, label in zip(axes.flat, labels):
        ax.text(
        0.05, 0.95, label,
        transform=ax.transAxes,
        fontsize=14,
        #fontweight="bold",
        va="top",
        ha="left"
        )

#axes[0, 0].set_title(
#    r"$B_z = 0.001\ \mathrm{a.u.}$",
#    fontsize=18,
#    pad=8
#)
#axes[0, 1].set_title(
#    r"$B_z = 0.01\ \mathrm{a.u.}$",
#    fontsize=18,
#    pad=8
#)
#axes[0, 2].set_title(
#    r"$B_z = 0.1\ \mathrm{a.u.}$",
#    fontsize=15,
#    pad=8
#)
for ax in axes.flat:
    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", labelsize=14)
for ax in axes.flat:
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#plt.legend()
plt.tight_layout()
plt.savefig("GOE.png", dpi=300, bbox_inches="tight")




        
