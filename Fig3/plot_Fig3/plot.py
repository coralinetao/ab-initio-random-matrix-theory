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

labels = ['(a)', '(b)', '(c)', '(d)', '(e)',
 '(f)', '(g)', '(h)', '(i)', '(j)']

plot_files = [
    # first row: no "noee"
"ch3_oxi_adz_3000.txt","alanine_adz_4000.txt",
"ph_nh2_ch3_adz.txt","benzene1_adz_B0.txt",
"c6h6_adz_4000.txt",
    # second row: with "noee"
    "ch3_oxi_acqz_noee.txt",
    "alanine_acqz_noee.txt",
    "ph_nh2_ch3_acqz_noee.txt",
    "benzene_perturbed_acqz_noee.txt",
    "benzene_eq_acqz_noee.txt",
]

#colors = ["black","black","black","black","black","tab:blue","tab:blue","tab:blue","tab:blue","tab:blue"]
colors = ["tab:blue","tab:blue","tab:blue","tab:blue","tab:blue","tab:blue","tab:blue","tab:blue","tab:blue","tab:blue"]
bin_width2 = 0.1
s_plot = np.linspace(0, 10, 400)
P_GOE = wigner_dyson(s_plot, np.pi/4.0)
P_POS = poisson_DOS(s_plot)
fig, axes = plt.subplots(2, 5, figsize=(15, 5), sharex=True)


for i in range(2):
    for j in range(5):
        idx = i * 5 + j
        x_centers2, hist_y2 = np.loadtxt(
            plot_files[idx],
            delimiter=",",
            skiprows=1,
            unpack=True
        )

        ax = axes[i, j]
        ax.bar(x_centers2, hist_y2, width=bin_width2, align='center', color=colors[idx],edgecolor="k",linewidth=0.1)
        if i == 0 and j == 0:
            ax.plot(s_plot, P_GOE, color="r", label=r"Wigner–Dyson (GOE), $P(s)=\frac{\pi}{2}s\,e^{-\pi s^2/4}")
            ax.plot(s_plot, P_POS, color="orange", label=r"Poisson, $P(s)=e^{-s}$")
        else:
            ax.plot(s_plot, P_GOE, color="r")
            ax.plot(s_plot, P_POS, color="orange")

        ax.set_xlim(-0.2, 6)
        #ax.set_ylim(0,max(hist_y2)+0.1)
        ax.text(
            0.85, 0.9,
            labels[idx],
            transform=ax.transAxes,
            fontsize=14,
       #     fontweight="bold"
        )



#axes[1,1].set_xlabel("Unfolded level spacing s = ΔE / ⟨ΔE⟩")
#axes[1,0].set_xlabel("Unfolded level spacing s = ΔE / ⟨ΔE⟩")
#axes[0,0].set_ylabel("Probability density P(s)")
#axes[1,0].set_ylabel("Probability density P(s)")

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

#plt.tight_layout(rect=[0, 0, 1, 0.92])  # leave room for legend
#plt.show()

for ax in axes.flat:
    ax.tick_params(axis="x", labelsize=12)
    ax.tick_params(axis="y", labelsize=12)
for ax in axes.flat:
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#plt.show()
plt.tight_layout()

fig.subplots_adjust(hspace=0.1)
plt.savefig("GOE.png", dpi=300, bbox_inches="tight")




        
