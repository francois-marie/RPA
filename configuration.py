import numpy as np
from numpy import pi
import matplotlib.pyplot as plt

class Configuration():

    def __init__(self, n_k, n_q, t, beta, mu, U):
        self.n_k = n_k
        self.t = t
        self.beta = beta
        self.mu = mu
        self.U = U
        k_lin_mesh = np.linspace(-pi, pi, n_k, endpoint=False)
        self.kx, self.ky = np.meshgrid(k_lin_mesh, k_lin_mesh)
        q_lin_mesh = np.linspace(-pi, pi, n_q, endpoint=True)
        self.qx, self.qy = np.meshgrid(q_lin_mesh, q_lin_mesh)
        self.energy = self._get_energy()
        self.fermi = self._get_fermi()
        self.chi0 = np.vectorize(self._get_chi0)(self.qx, self.qy)
        self.chiRPA = self.chi0/(1-self.U*self.chi0)

    def _get_energy(self):
        return (epsilon(self.kx, self.ky, self.t, self.mu))

    def _get_fermi(self):
        return (fermi(self.energy, self.beta))


    def chi0_term(self, qx, qy, kx, ky, t=1., mu=0., beta=15.):
        eps_kq, eps_k = epsilon(kx + qx, ky + qy, t, mu), epsilon(kx, ky, t, mu)

        mask = abs(eps_k - eps_kq) > 1e-12
        nmask = np.invert(mask)
        A = np.zeros_like(kx)
        A[mask] = (-(fermi(eps_k[mask], beta) - fermi(eps_kq[mask], beta)) / (eps_k[mask] - eps_kq[mask]))
        A[nmask] = der_fermi(eps_k[nmask], beta)
        return (A)


    def _get_chi0(self, qx, qy):
        return 1. / self.n_k**2 * np.sum(self.chi0_term(qx, qy, self.kx, self.ky, self.t, self.mu, self.beta))


def config_to_image(config, grandeur, xlabel, ylabel, title):
    """Turn an array into an image"""
    zmin, zmax = grandeur.min(), grandeur.max()
    plt.contourf(config.kx, config.ky, grandeur, cmap="terrain", levels=np.linspace(zmin, zmax, 100))
    plt.colorbar(ticks=np.asarray(np.linspace(zmin, zmax, 6)).round(decimals=2)).ax.tick_params(labelsize=15)
    plt.xlim(-pi, pi); plt.ylim(-pi, pi)
    plt.xlabel(r'$'+xlabel+'$', fontsize=15); plt.ylabel(r'$'+ylabel+'$', fontsize=15)
    plt.title(r"$"+title+"$", fontsize=25)
    plt.axis('equal')
    plt.show()

def epsilon(kx, ky, t=1., mu=0.):
    return -2 * t * (np.cos(kx) + np.cos(ky)) - mu

def fermi(e, beta=15.):
    return 1. / (np.exp(beta*e) + 1.)

def der_fermi(e, beta=15.):
    return beta * np.exp(beta * e)*fermi(e, beta)**2

if __name__ == "__main__":

    config = Configuration(40, 40, 1., 15., 0., 2.)
    config_to_image(config, config.chiRPA, 'q_x', 'q_y', '\chi_{RPA}(\mathbf{q},\Omega=0)')
