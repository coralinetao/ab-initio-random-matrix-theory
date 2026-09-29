# Figure 3 – Many-Particle Level-Spacing Statistics

This directory contains the data and plotting scripts used to generate
Figure 3.

Figure 3 compares nearest-neighbor level-spacing statistics for
interacting many-particle electronic excitations with a corresponding
noninteracting reference spectrum.

## Interacting spectrum

The interacting many-particle spectra are obtained from singlet
configuration-interaction singles (CIS) calculations using the
aug-cc-pVDZ basis set.

For each molecular system, a selected energy window of the CIS
spectrum is unfolded and used to calculate the nearest-neighbor
level-spacing distribution.

The resulting distributions are compared with the random-matrix
Wigner-Dyson distribution and the Poisson distribution.

## Noninteracting spectrum

The noninteracting reference spectra are constructed from the
Hartree-Fock orbital energies as

$$
E_{ia} = E_{\mathrm{gs}} + \epsilon_a - \epsilon_i,
$$


where \(i\) and \(a\) denote occupied and virtual Hartree-Fock
orbitals, respectively.

These spectra therefore contain the same underlying single-particle
orbital structure but exclude the interaction-induced mixing present
in the CIS calculation.

The construction, unfolding, and level-spacing analysis of the
noninteracting spectra are performed by `WD_HF.py` in the Figure 1
directory. The resulting `_noee` files are used here to generate the
noninteracting panels of Figure 3.

## Spectral unfolding and level statistics

For both interacting and noninteracting spectra, the selected energy
levels are unfolded using a smooth polynomial fit to the integrated
density of states.

Nearest-neighbor spacings are calculated from the unfolded energies
and normalized to unit mean spacing,

\[
s_n = \epsilon_{n+1}^{\mathrm{unfolded}}
      - \epsilon_n^{\mathrm{unfolded}}.
\]

The resulting spacing distributions are compared with the Poisson and
Wigner-Dyson reference distributions.

## Figure organization

- Figure 3(a–e): interacting singlet CIS spectra
- Figure 3(f–j): noninteracting many-particle reference spectra

## Reproduction

The plotting script in this directory reads the processed
level-spacing data for the interacting and noninteracting spectra and
generates the panels of Figure 3.

The noninteracting `_noee` data are generated from the Hartree-Fock
orbital energies using the analysis workflow provided with Figure 1.
