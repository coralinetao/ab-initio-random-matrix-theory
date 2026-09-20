# Figure 7 – Parametric Level Statistics

This directory contains the data and scripts used to generate Figure 7.

The workflow has three stages:

1. `Ex/`, `Ey/`, and `Ez/` independently calculate level velocities and
   curvatures from the field-dependent molecular electronic spectra.

2. `Fig7f/` performs the separate random-matrix calculation used for
   Figure 7(f).

3. `plot_Fig7/` contains copies of the processed outputs from these
   calculations together with `plot_Fig7.py`, which generates the final
   Figure 7 panels.
