from __future__ import absolute_import
from pytraj.datasets.DataSet_MatrixDbl import DataSet_MatrixDbl
from .base import plt
import numpy as np

def plot_matrix(dset, *args, **kwd):
    """plot matrix and return pyplot object
    Parameters
    ---------
    dset: {DataSet_MatrixDbl, ...}

    Return: a tuple of (Figure object, AxesSubplot object, AxesImage object)
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)

    if isinstance(dset, DataSet_MatrixDbl):
        # get matrix data
        # need to reshape since dset stores data in 1D
        mat = np.asarray(dset.get_full_matrix()).reshape(dset.n_rows, dset.n_cols)
    if isinstance(dset, np.ndarray) and len(dset.shape) == 2:
        mat = dset
    cax = ax.matshow(mat, *args, **kwd)
    return (fig, ax, cax)
