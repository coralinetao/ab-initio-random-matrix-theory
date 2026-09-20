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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import FormatStrFormatter

def wigner_dyson(s, B):
    return 2 * B * s * np.exp(-B * s**2)

def wigner_DOS(s):
    return 0.5 / np.pi*np.sqrt(4.0-s**2) 

def poisson_DOS(s):
    return np.exp(-s)


main_labels = ["(a)", "(b)", "(c)", "(d)"]
inset_labels = ["(e)", "(f)", "(g)", "(h)"]

small_files = [
    ["36_9N_dz_nocc.txt",    "36_v6_10_dz_nocc.txt"],
    ["18_7N_dz_nocc.txt",    "18_v6_5_dz_nocc.txt"]
]

large_files = [
    ["36_9N_dz_full.txt",   "36_v6_10_dz_full.txt"],
    ["18_7N_dz_full.txt",    "18_v6_5_dz_full.txt"]
]

colors = [
    ["black",   "black"],
    ["tab:blue",    "tab:blue"]
]

bin_width2 = 0.08
s_plot = np.linspace(0, 10, 400)
P_GOE = wigner_dyson(s_plot, np.pi/4.0)
P_POS = poisson_DOS(s_plot)
fig, axes = plt.subplots(2, 2, figsize=(10, 6), sharex=True, sharey=True)

count = 0
for i in range(2):
    for j in range(2):
        x_centers2, hist_y2 = np.loadtxt(
            small_files[i][j],
            delimiter=",",
            skiprows=1,
            unpack=True
        )

        ax = axes[i, j]
        ax.bar(x_centers2, hist_y2, width=bin_width2, align='center', color=colors[i][j],edgecolor="k",linewidth=0.2)
        if i == 0 and j == 0:
            ax.plot(s_plot, P_GOE, color="r", label="Wigner–Dyson (GOE)")
            ax.plot(s_plot, P_POS, color="orange", label="Poisson")
        else:
            ax.plot(s_plot, P_GOE, color="r")
            ax.plot(s_plot, P_POS, color="orange")

        ax.set_xlim(-0.2, 10)

        ax.text(
            0.03, 0.9,
            main_labels[count],
            transform=ax.transAxes,
            fontsize=14,
       #     fontweight="bold"
        )

        x_in, y_in = np.loadtxt(
            large_files[i][j],
            delimiter=",",
            skiprows=1,
            unpack=True
        )

        inset = inset_axes(ax, width="50%", height="50%", loc="upper right")
        inset.bar(x_in, y_in, width=bin_width2, align='center', color=colors[i][j],edgecolor="k",linewidth=0.2)
        inset.plot(s_plot, P_GOE, color='r', linewidth=1)
        inset.plot(s_plot,P_POS,color='orange')
        inset.grid(False)
        # make inset cleaner
        inset.set_xlim(-0.2, 5)          # or tighter if you want a zoom
        inset.tick_params(labelsize=6)

        inset.text(
            0.05, 0.85,
            inset_labels[count],
            transform=inset.transAxes,
            fontsize=12,
       #     fontweight="bold"
        )

        count += 1


#axes[1,1].set_xlabel("Unfolded level spacing s = ΔE / ⟨ΔE⟩", fontsize=15)
#axes[1,0].set_xlabel("Unfolded level spacing s = ΔE / ⟨ΔE⟩", fontsize=15)
#axes[0,0].set_ylabel("Probability density P(s)", fontsize=15)
#axes[1,0].set_ylabel("Probability density P(s)", fontsize=15)

fig.supxlabel(r"Unfolded level spacing $s$" , fontsize=18)
fig.supylabel(r"Probability density $P(s)$", fontsize=18)


#handles, labels = axes[0,0].get_legend_handles_labels()
#
#fig.legend(
#    handles, labels,
#    loc="upper center",
#    ncol=2,          # put GOE and Poisson side by side
#    frameon=False
#)

#plt.tight_layout(rect=[0, 0, 1, 0.92])  # leave room for legend
#plt.show()
    
for ax in axes.flat:
    ax.tick_params(axis="x", labelsize=12)
    ax.tick_params(axis="y", labelsize=12)
for ax in axes.flat:
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.show()
plt.tight_layout()
#plt.savefig("helicene.png", dpi=300, bbox_inches="tight")




        
