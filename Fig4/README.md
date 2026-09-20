# Figure 4 – Spectral Form Factor

This directory contains the unfolded many-particle energy spectra and
analysis script used to generate the spectral form factor results in
Figure 4.

## Spectral form factor

The spectral form factor (SFF) is calculated from the unfolded
many-particle energy levels as

\[
K(\tau)
=
\frac{1}{M^2}
\left|
\sum_{n=1}^{M}
e^{-2\pi i \tau \epsilon_n}
\right|^2,
\]

where \(\epsilon_n\) are the unfolded energy levels within a spectral
window containing \(M\) states.

The dimensionless time is

\[
\tau = \frac{t}{t_H},
\]

where \(t_H\) is the Heisenberg time. Because the spectra are unfolded
to unit mean level spacing, \(\Delta=1\), the Heisenberg time is

\[
t_H = \frac{2\pi}{\Delta} = 2\pi
\]

in atomic units. This gives the factor \(2\pi\) in the Fourier phase
used in `plot_SFF.py`.

With the normalization used here,

\[
K(0)=1,
\]

and the long-time plateau for a window containing \(M\) levels is
\(1/M\).

## Window averaging

To reduce fluctuations in the spectral form factor, the calculation is
performed over overlapping windows of the unfolded spectrum.

The parameters used for Figure 4 are

- window size: 300 consecutive energy levels
- window step: 50 energy levels

For each window, `plot_SFF.py` calculates \(K(\tau)\). The final
spectral form factor is obtained by averaging over all windows within
the selected unfolded spectrum.

## Input files

The `.dat` files contain the unfolded many-particle energy levels for
the molecular systems analyzed in Figure 4:

- alanine
- symmetry-broken \(C_1\) benzene
- \(D_{6h}\) benzene
- methyloxirane
- 1-phenylethylamine

These unfolded spectra are obtained from the same singlet CIS
electronic excitation energies used for the many-particle
level-statistics analysis.

## `plot_SFF.py`

`plot_SFF.py`:

1. reads and sorts each unfolded spectrum,
2. divides the spectrum into overlapping 300-state windows,
3. calculates the spectral form factor for each window,
4. averages the spectral form factor over all windows, and
5. plots the resulting \(K(\tau)\) curves on a logarithmic scale.

The horizontal dashed line at

\[
K = \frac{1}{300}
\]

indicates the long-time plateau expected from the normalization used
for a 300-state window.

## Running the script

Place the unfolded `.dat` files in the same directory as
`plot_SFF.py` and run

```bash
python plot_SFF.py
