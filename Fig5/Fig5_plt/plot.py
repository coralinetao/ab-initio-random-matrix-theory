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


files = [
    ["alanine_adz_B0.001.txt", "alanine_adz_B0.01.txt", "alanine_adz_B0.1.txt"],
#    ["c4h4_Bz0.001_adz.txt",    "c4h4_Bz0.01_adz.txt",    "c4h4_Bz0.1_adz.txt"],
    ["benzene_adz_Bz0.001.txt",
     "benzene_adz_Bz0.01.txt",
     "benzene_adz_Bz0.1.txt"],
]
bin_width2 = 0.05
s_plot = np.linspace(0, 10, 400)
P_GOE = wigner_dyson(s_plot, np.pi/4.0)
P_GUE = wigner_dyson_GUE(s_plot)

fig, axes = plt.subplots(2, 3, figsize=(10, 5), sharex=True, sharey=True)

for i in range(2):
    for j in range(3):
        x_centers2, hist_y2 = np.loadtxt(
            files[i][j],
            delimiter=",",
            skiprows=1,
            unpack=True
        )

        ax = axes[i, j]
        ax.bar(x_centers2, hist_y2, width=bin_width2, align='center', edgecolor='k',color="tab:blue",linewidth=0.1)
        ax.plot(s_plot,P_GOE,color='r',label="Wigner–Dyson (GOE)")
        ax.plot(s_plot,P_GUE,color='purple' ,label="Wigner–Dyson (GUE)")
        ax.set_xlim(-0.2, 5)

#        if i == 0 and j == 0:
#            ax.plot(s_plot, P_GOE, color="r", label=r"Wigner–Dyson (GOE), $P(s)=\frac{\pi}{2}s\,e^{-\pi s^2/4}$")
#            ax.plot(s_plot,P_GUE,color='purple' ,label=r"Wigner–Dyson (GUE), $P_{\mathrm{GUE}}(s)=\frac{32}{\pi^2}s^2\,e^{-\frac{4}{\pi}s^2}$")
#        else:
#            ax.plot(s_plot, P_GOE, color="r")
#            ax.plot(s_plot,P_GUE,color='purple')

fig.supxlabel(r"Unfolded level spacing $s$", fontsize=15)
fig.supylabel(r"Probability density $P(s)$", fontsize=15)



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
        0.02, 0.95, label,
        transform=ax.transAxes,
        fontsize=15,
        #fontweight="bold",
        va="top",
        ha="left"
        )

axes[0, 0].set_title(
    r"$B_z = 0.001\ \mathrm{a.u.}$",
    fontsize=15,
    pad=8
)
axes[0, 1].set_title(
    r"$B_z = 0.01\ \mathrm{a.u.}$",
    fontsize=15,
    pad=8
)
axes[0, 2].set_title(
    r"$B_z = 0.1\ \mathrm{a.u.}$",
    fontsize=15,
    pad=8
)

for ax in axes.flat:
    ax.tick_params(axis="x", labelsize=11)
    ax.tick_params(axis="y", labelsize=11)
for ax in axes.flat:
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#plt.legend()
plt.tight_layout()
#plt.savefig("bfield_final.png", dpi=300, bbox_inches="tight")
plt.show()




        
