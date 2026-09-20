# Figure 1 – Single-Particle Level-Spacing Statistics

This directory contains the scripts and data used for the
single-particle level-spacing analysis in Figure 1.

The same Hartree-Fock orbital energies are also used to construct the
non-interacting many-particle reference spectra shown in Figure 3 and
the density-of-states analysis shown in Supplemental Figure S1.

Each subdirectory contains a Python script together with the input data
required to run that script.

## WD_HF

`WD_HF.py` performs the statistical analysis starting from the
Hartree-Fock orbital energies.

For the single-particle spectrum, the script:

- selects the orbital-energy window used for the analysis,
- performs polynomial spectral unfolding,
- calculates nearest-neighbor level-spacing and gap-ratio statistics,
- generates the processed spacing data used for Figure 1, and
- generates the density-of-states data used for Supplemental Figure S1.

The script also constructs the non-interacting many-particle
single-excitation spectrum from the Hartree-Fock orbital energies as

E_ia = E_gs + epsilon_a - epsilon_i,

where i and a denote occupied and virtual Hartree-Fock orbitals,
respectively.

The non-interacting spectrum is unfolded and analyzed using the same
level-statistics procedure. The resulting processed data are used for
the non-interacting results in Figure 3.

## plot_Fig1

Contains the processed single-particle level-spacing data and plotting
script used to generate Figure 1.

## plot_FigS1

Contains the processed single-particle density-of-states data and
plotting script used to generate Supplemental Figure S1.
