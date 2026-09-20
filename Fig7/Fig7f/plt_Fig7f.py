import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# Random-matrix ensembles
# ============================================================
def goe(N, rng):
    """
    Generate an N x N real-symmetric random matrix (GOE).
    """
    A = rng.normal(size=(N, N))
    return (A + A.T) / 2.0

def gue(N, rng):
    """
    Generate an N x N complex-Hermitian random matrix (GUE).
    """
    A = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))
    return (A + A.conj().T) / 2.0

def H_crossover(H0, H1, B):
    """
    Construct a GOE-to-GUE crossover Hamiltonian

        H(B) = (H0 + B H1) / sqrt(1 + B^2).

    At B = 0 the Hamiltonian is GOE. Increasing |B| introduces
    the complex-Hermitian contribution H1 and breaks
    time-reversal symmetry.
    """
    return (H0 + B * H1) / np.sqrt(1.0 + B**2)

# ============================================================
# Spectral quantities
# ============================================================
def bulk_indices(N, frac=0.6):
    """
    Return indices corresponding to the central fraction of the
    spectrum. Edge states are excluded to reduce finite-size and
    density-of-states effects.
    """
    m = int((1 - frac) * N / 2)
    return np.arange(m, N - m)

def mean_spacing(evals, idx):
    """
    Mean nearest-neighbor level spacing in the selected bulk window.
    """
    return np.mean(np.diff(evals[idx]))

# ============================================================
# Curvature statistic
# ============================================================
def curvature(E_minus, E0, E_plus, h):
    return (E_plus - 2.0 * E0 + E_minus) / h**2

    """
    Compute the ensemble-averaged squared level curvature

        <(K / Delta)^2>

    as a function of the GOE-GUE crossover parameter B.

    K is the second derivative of each energy level with respect
    to B, and Delta is the mean level spacing in the bulk of the
    spectrum.
    """
def compute_K2(N, nreal, B_positive, h, seed=0):
    rng = np.random.default_rng(seed)
    idx = bulk_indices(N)

    K2 = np.zeros(len(B_positive))

    # Average over independent random-matrix realizations
    for r in range(nreal):
        H0 = goe(N, rng)
        H1 = gue(N, rng)

        for i, B in enumerate(B_positive):

            # Eigenvalues at B-h, B, and B+h are used to evaluate
            # the level curvature by finite differences.
            Em = np.linalg.eigvalsh(H_crossover(H0, H1, B - h))
            E0 = np.linalg.eigvalsh(H_crossover(H0, H1, B))
            Ep = np.linalg.eigvalsh(H_crossover(H0, H1, B + h))

            K = curvature(Em[idx], E0[idx], Ep[idx], h)

            # Normalize the curvature by the local mean
            # level spacing.
            Delta = mean_spacing(E0, idx)
            K /= Delta

            K2[i] += np.mean(K**2)
    # Ensemble average
    K2 /= nreal
    return K2

# ============================================================
# Main calculation
# ============================================================
if __name__ == "__main__":

    # Random-matrix dimension
    N = 250

    # Number of independent random-matrix realizations
    nreal = 80

    # Finite-difference step used to calculate curvature
    h = 2e-2

    # Maximum GOE-GUE crossover parameter
    Bmax = 4.0

    # Because the statistic is symmetric under B -> -B,
    # calculate only the positive-B branch.
    B_pos = np.linspace(0, Bmax, 30)
    K2_pos = compute_K2(N, nreal, B_pos, h)

    # Reflect the positive-B result to obtain the full
    # symmetric curve.
    B_full = np.concatenate((-B_pos[::-1], B_pos))
    K2_full = np.concatenate((K2_pos[::-1], K2_pos))

    # ========================================================
    # Save processed data
    # ========================================================
    data = np.column_stack((B_full, K2_full))
    np.savetxt(
        "curvature_vs_B.txt",
        data,
        header="B    <(K/Delta)^2>",
        comments=''
    )

    # ========================================================
    # Plot
    # ========================================================
    plt.figure()
    plt.plot(B_full, K2_full)   # no markers
    plt.xlabel("B")
    plt.ylabel(r"$\langle (K/\Delta)^2 \rangle$")
    plt.title("GOE → GUE curvature crossover")
    plt.tight_layout()
    plt.show()
