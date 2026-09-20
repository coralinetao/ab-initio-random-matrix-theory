import ast
import numpy as np
import matplotlib.pyplot as plt

with open("ravg.txt", "r") as f:
    ravg1 = np.array(ast.literal_eval(f.readline()), dtype=float)
    ravg2 = np.array(ast.literal_eval(f.readline()), dtype=float)

steps = np.arange(1, len(ravg1) + 1)

fig, ax = plt.subplots(figsize=(3.4, 2.8))

ax.plot(
    steps,
    ravg1,
    marker="o",
    markersize=5,
    linewidth=1.8,
)

ax.plot(
    steps,
    ravg2,
    marker="s",
    markersize=5,
    linewidth=1.8)

ax.axhline(
    0.3863,
    linestyle="--",
    linewidth=1.2,
    color = 'yellow'
)

ax.axhline(
    0.5359,
    linestyle="--",
    linewidth=1.2,
    color = 'red'
)

ax.set_xlabel(
    r"Interpolated geometry ($0.0138\ \mathrm{\AA}$ per step)",
    fontsize=8
)

ax.set_ylabel(
    r"$\langle r\rangle$",
    fontsize=10
)

ax.set_xticks(steps)

ax.tick_params(
    axis="both",
    which="major",
    labelsize=11,
    direction="in",
    length=5,
    width=1.0,
    top=True,
    right=True
)

ax.set_ylim(0.1, 0.55)

ax.legend(
    fontsize=9,
    frameon=False
)

fig.tight_layout()

#fig.savefig(
#    "ravg_vs_step.pdf",
#    bbox_inches="tight"
#)
#
#fig.savefig(
#    "ravg_vs_step.png",
#    dpi=600,
#    bbox_inches="tight"
#)

plt.show()
