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
"ch3_oxi_adz_3000_DOS.txt","alanine_adz_4000_DOS.txt",
"ph_nh2_ch3_adz_DOS.txt","benzene1_adz_B0_DOS.txt",
    # second row: with "noee"
    "ch3_oxi_acqz_noee_DOS.txt",
    "alanine_acqz_noee_DOS.txt",
    "ph_nh2_ch3_acqz_noee_DOS.txt",
    "benzene_perturbed_acqz_noee_DOS.txt",
]

#colors = ["black","black","black","black","black","tab:blue","tab:blue","tab:blue","tab:blue","tab:blue"]
colors = ["tab:blue","tab:blue","tab:blue","tab:blue","tab:blue","tab:blue","tab:blue","tab:blue","tab:blue","tab:blue"]
bin_width2 = 0.1
s_plot = np.linspace(0, 10, 400)
P_GOE = wigner_dyson(s_plot, np.pi/4.0)
P_POS = poisson_DOS(s_plot)
fig, axes = plt.subplots(2, 4, figsize=(15, 7), sharey="row")


for i in range(2):
    for j in range(4):
        idx = i * 4 + j
        x_centers2, hist_y2 = np.loadtxt(
            plot_files[idx],
            delimiter=",",
            skiprows=1,
            unpack=True
        )
        with open(plot_files[idx], "r") as f:
             header_line = f.readline().strip()
        val0_str, val1_str = header_line.split("_")
        cut0 = float(val0_str)
        cut1 = float(val1_str)


        ax = axes[i, j]
        ax.bar(x_centers2, hist_y2, width=bin_width2, align='center', color=colors[idx],edgecolor="k",linewidth=0.01)
        #ax.set_ylim(0,max(hist_y2)+0.1)
        ax.axvspan(cut0, cut1, color='orange', alpha=0.2)
        if i > 0:
           nonzero_indices = np.where(hist_y2 > 0)[0]
           ax.set_xlim(x_centers2[nonzero_indices[0]]-10,x_centers2[nonzero_indices[0]]+100)
           ax.set_ylim(0,0.2)
        ax.text(
            0.85, 0.9,
                   labels[idx],
            transform=ax.transAxes,
            fontsize=16,
       #     fontweight="bold"
        )



#axes[1,1].set_xlabel("Unfolded level spacing s = ΔE / ⟨ΔE⟩")
#axes[1,0].set_xlabel("Unfolded level spacing s = ΔE / ⟨ΔE⟩")
#axes[0,0].set_ylabel("Probability density P(s)")
#axes[1,0].set_ylabel("Probability density P(s)")

fig.supxlabel(r"Excited state energies $E$ (Hartree)", fontsize=18)
fig.supylabel(r"Probability density $P(E)$", fontsize=18)

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
    ax.tick_params(axis="x", labelsize=13)
    ax.tick_params(axis="y", labelsize=13)
for ax in axes.flat:
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
#plt.show()
plt.tight_layout()

fig.subplots_adjust(hspace=0.1)
plt.savefig("MP_DOS.png", dpi=300, bbox_inches="tight")




        
