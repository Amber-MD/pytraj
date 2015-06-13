from __future__ import absolute_import
from pytraj.datasets.DatasetMatrixDouble import DatasetMatrixDouble
from ..utils import _import_numpy
from .base import plt, np

def plot_matrix(dset, *args, **kwd):
    """plot matrix and return pyplot object
    Parameters
    ---------
    dset: {DataSet_MatrixDbl, ...}

    Return: a tuple of (Figure object, AxesSubplot object, AxesImage object)
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)

    #if isinstance(dset, DataSet_MatrixDbl):
    # use below for all kind of matrix datatypes
    if 'matrix' in dset.dtype.lower(): 
        # get matrix data
        # need to reshape since dset stores data in 1D
        mat = np.asarray(dset.get_full_matrix()).reshape(dset.n_rows, dset.n_cols)
    if isinstance(dset, np.ndarray) and len(dset.shape) == 2:
        mat = dset
    cax = ax.matshow(mat, *args, **kwd)
    return (fig, ax, cax)
