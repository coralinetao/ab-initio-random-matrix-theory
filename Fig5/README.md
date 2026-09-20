# Figure 5 – Magnetic-Field-Induced GOE-to-GUE Crossover

This directory contains the electronic excitation energies, analysis
scripts, and processed data used to generate Figure 5.

Figure 5 examines the change in many-particle level statistics when
time-reversal symmetry is broken by an external magnetic field.

## Electronic-structure calculations

The many-particle electronic excitation energies are obtained from
singlet configuration-interaction singles (CIS) calculations using the
aug-cc-pVDZ basis set.

The zero-field spectra exhibit level statistics associated with the
Gaussian Orthogonal Ensemble (GOE) when the relevant spatial
symmetries are broken.

An external magnetic field breaks time-reversal symmetry, producing a
crossover toward the Gaussian Unitary Ensemble (GUE).

## Spectral unfolding

For each magnetic-field strength, a selected portion of the
many-particle spectrum is unfolded using a smooth polynomial fit to
the integrated density of states.

If the ordered energy levels are

\[
E_1 < E_2 < \cdots,
\]

the spectral staircase is defined by

\[
N(E_i)=i.
\]

A smooth approximation \(\overline{N}(E)\) is fitted to the staircase,
and the unfolded energies are

\[
\epsilon_i=\overline{N}(E_i).
\]

The nearest-neighbor spacings are then

\[
s_i=\epsilon_{i+1}-\epsilon_i.
\]

The unfolding procedure is the same as that used for the
many-particle level-spacing analysis in Figure 3.

## Level-spacing statistics

The nearest-neighbor spacing distributions are compared with the
random-matrix predictions for the Gaussian Orthogonal Ensemble (GOE)
and Gaussian Unitary Ensemble (GUE).

For time-reversal-invariant systems, the relevant Wigner-Dyson
reference distribution is the GOE form,

\[
P_{\mathrm{GOE}}(s)
=
\frac{\pi}{2}s
\exp\left(-\frac{\pi s^2}{4}\right).
\]

When time-reversal symmetry is broken by the magnetic field, the
corresponding GUE reference distribution is

\[
P_{\mathrm{GUE}}(s)
=
\frac{32}{\pi^2}s^2
\exp\left(-\frac{4s^2}{\pi}\right).
\]

Figure 5 uses these distributions to visualize the change in spectral
statistics as the magnetic field is introduced.

## Input and processed data

The input `.dat` files contain the many-particle electronic excitation
energies obtained from the CIS calculations at the magnetic-field
strengths considered in Figure 5.

The accompanying `.param` files specify the parameters used for
selecting and unfolding the relevant portion of each spectrum.

The analysis scripts generate processed level-spacing distributions
that are subsequently used by the Figure 5 plotting script.

## Reproduction

Run the level-statistics analysis for the desired magnetic-field
calculation using its corresponding `.dat` and `.param` files.

The resulting processed spacing distributions can then be used by the
Figure 5 plotting script to reproduce the published figure.
