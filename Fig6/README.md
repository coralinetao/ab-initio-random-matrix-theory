# Figure 6 – Parametric Evolution of Electronic Energy Levels

This directory contains the data and scripts used to generate the two
level-evolution ("noodle") plots shown in Figure 6.

The two panels were generated separately and are therefore organized
into separate subdirectories:

- `Fig6a/` – data and script for Figure 6(a)
- `Fig6b/` – data and script for Figure 6(b)

The electronic energy levels are obtained from singlet
configuration-interaction singles (CIS) calculations using the
aug-cc-pVDZ basis set.

## Figure 6(a)

`Fig6a/` contains the data and plotting script used to generate the
level trajectories shown in Figure 6(a).

The script reads the parameter-dependent electronic energy levels and
plots the selected states as a function of the corresponding external
parameter.

## Figure 6(b)

`Fig6b/` contains the data and analysis script used to generate the
level trajectories shown in Figure 6(b).

The input file `E.dat` contains the electronic energy levels as a
function of the perturbation parameter:

- rows correspond to electronic states,
- columns correspond to successive values of the perturbation
  parameter.

At each parameter value, the spectrum is unfolded independently using
a polynomial fit to the smooth spectral staircase,

\[
N(E_i)=i,
\]

with unfolded energies defined as

\[
\epsilon_i=\overline{N}(E_i).
\]

A broad spectral window is used to determine the smooth density of
states, after which the subset of states displayed in Figure 6(b) is
selected.

### Rescaled parameter

The horizontal coordinate in Figure 6 is rescaled using the
characteristic level velocity of the unfolded spectrum.

The level velocities are calculated as

\[
V_i(t)=\frac{d\epsilon_i}{dt},
\]

and the velocity scale is defined by

\[
C(0)=\left\langle V_i(t)^2\right\rangle.
\]

The rescaled parameter is then

\[
x=\sqrt{C(0)}\,(t-t_0).
\]

The level velocities are therefore used only to determine the
horizontal scale of the plot; the quantities displayed in the figure
are the unfolded electronic energy levels \(\epsilon_i(x)\).

## Reproduction

Each panel can be reproduced independently by running the plotting
script in its corresponding directory:

```bash
cd Fig6a
python plot.py
