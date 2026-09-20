# Figure 2 – Real- and Momentum-Space Orbital Distributions

This directory contains the molecular-orbital cube files and plotting
script used to generate Figure 2.

## Files

`plot.py` reads a selected molecular-orbital cube file and generates
the corresponding real-space and momentum-space probability-density
plots.

`mo.<orbital>.cube` contains the molecular-orbital wavefunction on a
three-dimensional real-space grid.

The orbital to be analyzed is selected by setting `orb_num` in
`plot.py`.

## Real-space distribution

The real-space probability density is calculated as

\[
\rho(\mathbf{r}) = |\psi(\mathbf{r})|^2.
\]

The three-dimensional density is projected onto the \(x\)-\(y\) plane
by integrating over the \(z\) direction,

\[
\rho(x,y) = \int dz\,|\psi(x,y,z)|^2.
\]

This projected density is used for the real-space panels of Figure 2.

## Momentum-space distribution

The molecular orbital is transformed to momentum space using a
three-dimensional fast Fourier transform (FFT). The momentum-space
probability density is

\[
\rho(\mathbf{p}) = |\tilde{\psi}(\mathbf{p})|^2.
\]

For visualization, the three-dimensional momentum-space density is
projected onto the \(p_x\)-\(p_y\) plane by integrating over \(p_z\),

\[
\rho(p_x,p_y)
=
\int dp_z\,
|\tilde{\psi}(p_x,p_y,p_z)|^2.
\]

This projected distribution is used for the momentum-space panels of
Figure 2.

## Equimomentum contour

The red circle in each momentum-space panel is a reference
equimomentum contour defined by the average kinetic energy of the
orbital.

The script evaluates the average kinetic energy directly from the
three-dimensional momentum-space distribution,

\[
\langle T\rangle
=
\frac{\langle p^2\rangle}{2m_e},
\]

where

\[
\langle p^2\rangle
=
\langle p_x^2+p_y^2+p_z^2\rangle.
\]

For a one-electron orbital,

\[
\langle T\rangle
=
\epsilon-\langle \hat V_{\mathrm{eff}}\rangle.
\]

The radius of the reference contour is therefore determined by

\[
p_0^2
=
2m_e\langle T\rangle
=
2m_e\left(
\epsilon-\langle\hat V_{\mathrm{eff}}\rangle
\right).
\]

The circle displayed in the \(p_x\)-\(p_y\) plane satisfies

\[
p_x^2+p_y^2=p_0^2.
\]

Although \(p_z\) is not displayed explicitly, it is integrated over
in the plotted momentum-space density and is included in the full
three-dimensional \(\langle p^2\rangle\) used to determine the
reference radius.

## Running the script

Set the desired orbital number in `plot.py`, for example,

```python
orb_num = 861
