from __future__ import absolute_import
from .base import plt


def plot_matrix(dset, show=False, *args, **kwd):
    """plot matrix and return pyplot object
    Parameters
    ---------
    dset: {DataSet_MatrixDbl, ...}

    Return: a tuple of (Figure object, AxesSubplot object, AxesImage object)
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cax = ax.matshow(dset, *args, **kwd)
    if show:
        plt.show()
    return (fig, ax, cax)
