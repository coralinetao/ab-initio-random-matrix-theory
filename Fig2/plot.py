import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# ============================================================
# Momentum-space moments
# ============================================================

def momentum_moments_from_I3(I3, kx, ky, kz, m=1.0, hbar=1.0):
    """
    Calculate momentum moments and the average kinetic energy
    from the three-dimensional momentum-space probability density.

    Parameters
    ----------
    I3 : ndarray
        |psi(p)|^2 evaluated on a regular 3D momentum grid.
    kx, ky, kz : ndarray
        Wave-vector grids. Momentum is related to wave vector by
        p = hbar * k.
    m : float
        Electron mass. m = 1 in atomic units.
    hbar : float
        Reduced Planck constant. hbar = 1 in atomic units.

    Returns
    -------
    norm :
        Momentum-space normalization.
    p_mean :
        (<p_x>, <p_y>, <p_z>).
    p2_mean :
        <p^2> = <p_x^2 + p_y^2 + p_z^2>.
    p_rms :
        sqrt(<p^2>).
    ke_mean :
        Average kinetic energy <T> = <p^2> / (2m).
    """

    # Momentum-grid volume element.
    # The FFT grid is expressed initially in wave-vector coordinates k.
    dkx = float(kx[1] - kx[0])
    dky = float(ky[1] - ky[0])
    dkz = float(kz[1] - kz[0])
    d3k = dkx * dky * dkz

    # Construct the full 3D wave-vector grid.
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")

    # Convert wave vector to momentum: p = hbar * k.
    PX = hbar * KX
    PY = hbar * KY
    PZ = hbar * KZ

    # Momentum-space volume element:
    #
    #     d^3p = hbar^3 d^3k
    #
    d3p = (hbar**3) * d3k

    # Normalize the three-dimensional momentum-space density.
    norm = np.sum(I3) * d3p

    # First moments of the momentum distribution.
    # <p> = ∫ p |psi|^2 d^3p
    px_mean = np.sum(PX * I3) * d3p / norm
    py_mean = np.sum(PY * I3) * d3p / norm
    pz_mean = np.sum(PZ * I3) * d3p / norm

    # Total squared momentum:
    #
    #     p^2 = p_x^2 + p_y^2 + p_z^2
    #
    # Note that p_z is included here even though the final figure
    # displays a 2D projection onto the p_x-p_y plane.
    p2 = PX**2 + PY**2 + PZ**2
    p2_mean = np.sum(p2 * I3) * d3p / norm

    # Average kinetic energy of the orbital:
    #
    #     <T> = <p^2> / (2m)
    #
    ke_mean = p2_mean / (2.0 * m)

    return {
        "norm": norm,
        "p_mean": (px_mean, py_mean, pz_mean),
        "p2_mean": p2_mean,
        "p_rms": np.sqrt(p2_mean),
        "ke_mean": ke_mean,
    }


# ============================================================
# Read Gaussian cube file
# ============================================================
def read_cube(filename):
    """
    Read a molecular-orbital wavefunction from a Gaussian cube file.

    Returns the 3D orbital values, grid origin, grid dimensions,
    and the three real-space grid vectors.
    """
    with open(filename, "r") as f:
        title1 = f.readline()
        title2 = f.readline()

        # Number of atoms and origin of the real-space grid.
        parts = f.readline().split()
        natoms = int(parts[0])
        origin = np.array(parts[1:], float)

        # Grid dimensions and grid vectors along x, y, and z.
        nx, vx = parse_grid_line(f)
        ny, vy = parse_grid_line(f)
        nz, vz = parse_grid_line(f)
       
        # Skip atomic-coordinate lines.
        for _ in range(natoms):
            f.readline()

        # Read the molecular-orbital values on the 3D grid.
        data = []
        for line in f:
            data.extend([float(x) for x in line.split()])

    # Convert to 3D array
    cube = np.array(data).reshape((nx, ny, nz))
    return cube, origin, nx,ny,nz,vx, vy, vz


def parse_grid_line(f):
    """Read one grid-definition line from a cube file."""
    parts = f.readline().split()
    n = int(parts[0])
    vec = np.array(parts[1:], float)
    return n, vec


# ============================================================
# Select molecular orbital
# ============================================================
# Orbital number to analyze.
# The corresponding input file must be named:
#
#     mo.<orb_num>.cube
#
#orb_num = 704
#orb_num = 697
orb_num = 861 
phi, origin, nx, ny, nz, vx,vy,vz = read_cube("mo."+str(orb_num)+".cube")

# ============================================================
# Construct the real-space grid
# ============================================================

# Grid spacings along the Cartesian directions.
# The cube grids used here are aligned with x, y, and z.
dx = vx[0]
dy = vy[1]
dz = vz[2]
x = origin[0] + np.arange(nx)*dx
y = origin[1] + np.arange(ny)*dy
z = origin[2] + np.arange(nz)*dz

# Axis 2 corresponds to the z direction.
# Both real- and momentum-space plots below therefore show
# projections onto the x-y and p_x-p_y planes, respectively.
axis_plot = 2

print("dx: {}, dy: {}, dz: {} ".format(dx,dy,dz))

# ============================================================
# Real-space probability density
# ============================================================

# Three-dimensional real-space probability density:
#
#     rho(r) = |psi(r)|^2
#
Ix3 = np.abs(phi)**2  

# Check normalization of the orbital on the real-space grid.
dV = dx * dy * dz
norm_real = np.sum(Ix3) * dV
print("Real-space norm =", norm_real)

# Project the 3D real-space density onto the x-y plane by
# summing over the z grid.
#
# This represents the z-projected spatial distribution of the orbital.
Ixy3 = Ix3.sum(axis_plot) 

print("Ixy3 shape {}".format(Ixy3.shape))

# Plot the projected real-space probability density.
plt.figure()
plt.imshow(Ixy3.T, origin="lower", aspect="equal",
extent=[x[0], x[-1], y[0], y[-1]]#)
          , norm=LogNorm())
cbar=plt.colorbar(label=r"$\int|\psi(\mathbf{r})|^2 dr_z$")
cbar.ax.yaxis.label.set_size(18)
cbar.ax.tick_params(labelsize=15)
plt.xlabel(r"$x$ (bohr)", fontsize=18)
plt.ylabel(r"$y$ (bohr)", fontsize=18)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

plt.tight_layout()
plt.savefig(str(orb_num)+"_x.png", dpi=300, bbox_inches="tight")

# ============================================================
# Transform orbital from real space to momentum space
# ============================================================

# Three-dimensional Fourier transform of the real-space orbital.
#
# phi_k represents the orbital on the reciprocal-space grid.
phi_k = np.fft.fftn(phi, norm="ortho")          # 3D FFT
phi_k = np.fft.fftshift(phi_k)    # move k=0 to center
kx = np.fft.fftshift(np.fft.fftfreq(nx, d=dx)) * 2*np.pi
ky = np.fft.fftshift(np.fft.fftfreq(ny, d=dy)) * 2*np.pi
kz = np.fft.fftshift(np.fft.fftfreq(nz, d=dz)) * 2*np.pi

# ============================================================
# Momentum-space probability density
# ============================================================

# Three-dimensional momentum-space probability density.
I3 = np.abs(phi_k)**2   

# Calculate normalization, momentum moments, and average kinetic energy.
#
# Atomic units are used:
#
#     m_e = 1
#     hbar = 1
#
stats = momentum_moments_from_I3(I3, kx, ky, kz, m=1.0, hbar=1.0)

# Normalize the 3D momentum-space distribution.
I3/=stats["norm"]

# ============================================================
# Project momentum density onto the p_x-p_y plane
# ============================================================

# Sum over p_z to obtain the 2D distribution displayed in Fig. 2:
#
#     I(p_x,p_y) ~ integral dp_z |psi(p_x,p_y,p_z)|^2
#
# Thus p_z is not displayed explicitly in the figure. However,
# p_z IS included in the calculation of <p^2> and therefore in
# the average kinetic energy used to define the red reference circle.

I_xy = I3.sum(axis_plot)
print(stats)

# ============================================================
# Equimomentum reference contour
# ============================================================

# The red circle in Fig. 2 is an equimomentum reference contour
# in the displayed p_x-p_y plane.
#
# For a one-electron orbital,
#
#     <T> = epsilon - <V_eff>
#         = <p^2> / (2 m_e),
#
# where
#
#     <p^2> = <p_x^2 + p_y^2 + p_z^2>.
#
# The code evaluates <T> directly from the Fourier-transformed
# orbital rather than separately evaluating epsilon and <V_eff>.
#
# The radius of the reference contour is therefore defined by
#
#     p0^2 = 2 m_e <T>.
#
# In the displayed p_x-p_y plane, the circle satisfies
#
#     p_x^2 + p_y^2 = p0^2
#                   = 2 m_e (epsilon - <V_eff>).
#
# IMPORTANT:
# p_z has not been assumed to vanish. The plotted probability
# density is projected over p_z, while p_z contributes to the
# full 3D <p^2> used to determine the reference radius.

m = 1.0
hbar = 1.0
KE = stats["ke_mean"]

# k0 is the radius of the reference circle in reciprocal-space units.
k2 = 2.0*m*KE / hbar**2
k0 = np.sqrt(k2)
print(k2)

# ============================================================
# Momentum-space plots
# ============================================================

# Full projected momentum-space distribution.
plt.figure()
plt.imshow(
    I_xy.T,
    origin="lower",
    extent=[kx[0], kx[-1], ky[0], ky[-1]],
    aspect="equal",
)
plt.colorbar()

# Restrict the displayed momentum range for Fig. 2.
mask = (np.abs(kx) < 30)
kx_zoom = kx[mask]
ky_zoom = ky[mask]

I_zoom = I_xy[np.ix_(mask, mask)]

# Plot the momentum-space panel used in Fig. 2.
plt.figure()
plt.imshow(
    I_zoom.T,
    origin="lower",
    extent=[kx_zoom[0], kx_zoom[-1], ky_zoom[0], ky_zoom[-1]],
    aspect="equal"
)

# Draw the equimomentum reference circle:
#
#     p_x^2 + p_y^2 = 2 m_e <T>
#
theta = np.linspace(0, 2*np.pi, 800)
plt.plot(k0*np.cos(theta), k0*np.sin(theta), linewidth=1,color='red')
plt.xlabel(r"$p_x (\hbar/a_0)$", fontsize=18)
plt.ylabel(r"$p_y (\hbar/a_0)$", fontsize=18)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
cbar=plt.colorbar(label=r"$\int|\tilde\psi(\mathbf{p})|^2 dp_z$")
cbar.formatter.set_scientific(True)
cbar.formatter.set_powerlimits((0, 0))  # always use scientific notation
cbar.update_ticks()
cbar.ax.yaxis.label.set_size(18)
cbar.ax.tick_params(labelsize=15)
cbar.ax.yaxis.get_offset_text().set_size(15)
plt.tight_layout()
plt.savefig(str(orb_num)+"_p.png", dpi=300, bbox_inches="tight")
plt.show()

