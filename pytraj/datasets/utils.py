"""
Notes : might move to cython
"""
from __future__ import absolute_import
from ..DataSetList import DataSetList
from ..six_2 import zip

def stack(*args):
    """return a new DataSetList by joining (vstack)

    Parameters
    ----------
    *args : list/tuple of DataSetList

    Notes
    -----
        similiar to numpy.vstack

    Examples
    --------
        d1 = calc_dssp(traj1, dtype='dataset')
        d2 = calc_dssp(traj2, dtype='dataset')
        d3 = stack(d1, d2)
    """

    if len(args) <= 1:
        raise ValueError("need more than 1 DataSetList to stack")

    dslist0 = args[0].copy()
    for dslist in args[1:]:
        for d0, d in zip(dslist0, dslist):
            if d0.dtype != d.dtype:
                raise ValueError("Dont support stack different dtype together")
            d0.append(d)
    return dslist0
