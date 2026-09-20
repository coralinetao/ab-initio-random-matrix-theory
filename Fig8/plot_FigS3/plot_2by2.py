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

files = [["36_9N_dz_DOS.txt",
"36_v6_10_dz_DOS.txt"],["18_7N_dz_DOS.txt","18_v6_5_dz_DOS.txt"]]


#file_inset = "c6h6_adz_4000.txt"
bin_width1 = 0.1
bin_width2 = 0.01
s_plot = np.linspace(0, 10, 400)
P_GOE = wigner_dyson(s_plot, np.pi/4.0)
P_POS = poisson_DOS(s_plot)


fig, axes = plt.subplots(2, 2, figsize=(15, 8),constrained_layout=True )

for i in range(2):
    for j in range(2):
        x_centers2, hist_y2 = np.loadtxt(
            files[i][j],
            delimiter=",",
            skiprows=1,
            unpack=True
        )
        with open(files[i][j], "r") as f:
             header_line = f.readline().strip()
        val0_str, Nocc_str, val1_str = header_line.split("_")
        cut0 = float(val0_str)
        cut1 = float(val1_str)
        Nocc = float(Nocc_str)
        ax = axes[i, j]
        if i == 0 :
            ax.bar(x_centers2, hist_y2, width=bin_width1, align='center', edgecolor='k', linewidth=0.1,color='black')
        if i == 1 : 
            ax.bar(x_centers2, hist_y2, width=bin_width2, align='center', edgecolor='k', linewidth=0.01,color='tab:blue')
        #ax.set_xlim(-22.0, 100.0)
        #ax.set_ylim(0.0,1.2)
       # ax.axvline(x=cut0, color='r', linestyle='--')
        ax.axvline(x=Nocc, color='r', linestyle='--') 
        ax.axvspan(cut0, cut1, color='orange', alpha=0.2)

# Top row labels
#fig.text(0.5, 0.52, r"Orbital energies $\epsilon$ (Hartree)",
#         ha='center', va='center', fontsize=16)

fig.text(0.01, 0.75, r"Probability density $P(\epsilon)$",
         ha='center', va='center', rotation='vertical', fontsize=16)

# Bottom row labels
#fig.text(0.5, 0.01, r"Excited state energies $E$ (Hartree)",
#         ha='center', va='center', fontsize=16)

fig.text(0.01, 0.25, r"Probability density $P(E)$",
         ha='center', va='center', rotation='vertical', fontsize=16)

for j in range(2):
    axes[1, j].set_xlabel(r"Excited state energies $E$ (Hartree)", fontsize=16)

for j in range(2):
    axes[0, j].set_xlabel(r"Orbital energies $\epsilon$ (Hartree)", fontsize=16)

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
        fontsize=18,
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
    ax.tick_params(axis="x", labelsize=16)
    ax.tick_params(axis="y", labelsize=16)
for ax in axes.flat:
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#plt.legend()
ax.xaxis.set_major_formatter('{:.1f}'.format)


plt.tight_layout()
#plt.savefig("GOE.png", dpi=300, bbox_inches="tight")
plt.show()



        
