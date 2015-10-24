from matplotlib import pyplot as plt


def plot_phipsi(phi, psi, mark='ro'):
    """2D plot phi/psi
    Parameters
    ----------
    phi : ndarray, shape=(n, )

    psi : ndarray, shape=(n, )

    ptype : str, default='pmf' (2D free energy)
        plot type ('pmf', ...)
    """
    plt.plot(phi, psi, mark)
