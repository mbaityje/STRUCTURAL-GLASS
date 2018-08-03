"""
beylkin.py

Author: Ian S. Dunn, Columbia University Department of Chemistry
"""

# Import modules.
import sys
import numpy as np
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt
import scipy


class Beylkin:

    def __init__(self):

        self.x = None
        self.h = None
        self.N = None
        self.H = None
        self.ceigenvec = None
        self.gamma = None

    def load_data(self, x, h):

        assert(len(x) == len(h))

        self.N = int( (len(x) - 1) / 2 )
        self.x = x[:2*self.N+1]
        self.h = h[:2*self.N+1]

        self.T = self.x[-1] - self.x[0]

    def sample(self, f, start=0, end=1, N=214):

        # Number of points for sampling.
        self.N = N

        # Sampling grid.
        self.x = np.linspace(0, 1, 2 * N + 1, endpoint=True, dtype=np.complex128)

        # Evaluate function on grid.
        self.h = f(self.x)

    def plot_input(self):

        # Plot function.
        plt.plot(self.x, self.h)
        plt.show()

    def build_hankel(self):

        assert self.h is not None

        # Build Hankel matrix.
        self.H = np.zeros((self.N+1, self.N+1)).astype(np.complex128)
        for i in range(self.N + 1):
            for j in range(self.N + 1):
                self.H[i, j] = self.h[i + j]

        # Check that Hankel matrix is real and symmetric.
        assert np.allclose(self.H, self.H.T, atol=1.e-15)

    def eigen(self, nsing):

        assert self.H is not None

        # Initial Krylov vector.
        v0 = np.sin(np.arange(len(self.H)), dtype=np.complex128)

        # Calculate singular vectors and values of Hankel matrix.
        w, v = eigsh(self.H, k=nsing+1, v0=v0)

        self.ceigenvec = v[:, -1]

    def nodes(self):

        assert self.ceigenvec is not None

        # Evaluate roots of c-eigenpolynomial.
        gamma = np.roots(self.ceigenvec[::-1])

        # Remove large roots.
        max_gamma = (np.max(abs(self.h)) * 1000.) ** (1. / (2*self.N))
        large_inds = np.where(abs(gamma) > max_gamma)
        gamma = np.delete(gamma, large_inds)

        # Remove small roots.
        gamma = np.delete(gamma, np.where(abs(gamma) > 1.))

        # Store nodes.
        self.gamma = gamma

    def calculate_weights(self, nsing):

        assert self.gamma is not None
        assert self.gamma.size > 0

        # Build Vandermonde matrix from c-eigenroots.
        self.vand = np.vander(self.gamma, N=2*self.N+1).transpose()[::-1]

        # Normalize Vandermonde columns to improve conditioning of least squares by SVD.
        self.vand_norm = np.linalg.norm(self.vand, axis=0)
        self.vand /= self.vand_norm
        self.vand[np.where(abs(self.vand) < 1.e-14)] = 0

        # Change basis from complex exponentials to oscillating real exponentials.
        lamda_t = np.log(self.vand)
        omega_real_t = np.real(lamda_t)
        omega_imag_t = np.imag(lamda_t)
        self.vand = abs(self.vand)
        
        for i, wt in enumerate(omega_imag_t[1, :]):
        
            if wt > 1.e-14:
                self.vand[:, i] *= np.cos(omega_imag_t[:, i])
            elif wt < -1.e-14:
                self.vand[:, i] *= -np.sin(omega_imag_t[:, i])

        # Calculate Prony weights using least squares fit.
        lstsq_ret = scipy.linalg.lstsq(self.vand, self.h)
        self.weights = lstsq_ret[0]

        # Remove small weights.
        self.weights[np.where(abs(self.gamma) < 1.e-14)] = 0.

        # Sort weights.
        inds = np.argsort(abs(self.weights))

        # Set small weights to zero.
        self.weights[inds[:-nsing]] = 0

    def plot_components(self):

        # Plot absolute value of Prony components.
        for i in range(self.vand.shape[1]):
            plt.plot(abs(self.vand[:, i]))
        plt.show()

    def prony(self):

        assert self.vand is not None
        assert self.weights is not None

        # Construct approximate Prony approximation.
        self.approx = np.dot(self.vand, self.weights)

    def plot_prony(self):

        assert self.approx is not None

        # Plot original function and approximate Prony approximation.
        plt.plot(self.approx)
        plt.plot(self.h)
        plt.show()

    def correction(self):

        inds = np.where(abs(self.gamma) <= 1.)
        weights = np.copy(self.weights)
        #weights[inds] = 0.

        return np.dot(self.vand, weights)

    def driver(self, f, nsing):

        self.sample(f)

        self.build_hankel()

        self.eigen(nsing)

        self.nodes()

        self.calculate_weights(nsing)

        self.prony()

    def driver_load(self, x, h, nsing):

        self.load_data(x, h)

        self.build_hankel()

        self.eigen(nsing)

        self.nodes()

        self.calculate_weights(nsing)

        self.prony()

if __name__ == "__main__":

    f = bessel

    b = Beylkin()

    b.driver(f, int(sys.argv[1]))
