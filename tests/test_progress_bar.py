#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.testing import aa_eq
from pytraj.testing import tempfolder
import pytest

try:
    import tqdm
    from tqdm import tqdm_notebook
    import traitlets.traitlets.TraitError
except ImportError:
    tqdm = None

try:
    from pytraj.utils.progress import ProgressBarTrajectory
    import IPython
    has_ipython = True
except ImportError:
    IPython = None
    has_ipython = False

# local
from utils import fn


@unittest.skipUnless(has_ipython, 'need IPython')
def test_progress_log():
    """test_progress_log: simple test, just to make sure it is runnable
    """
    traj = pt.datafiles.load_tz2()[:10]

    p = ProgressBarTrajectory(traj, style='basic', every=5)
    pt.molsurf(p)

    p = ProgressBarTrajectory(traj, style='bar', every=5)
    pt.molsurf(p)

    p = ProgressBarTrajectory(traj, style='circle', every=5)
    pt.molsurf(p)

    # need to run on Jupyter notebook
    if tqdm is not None:
        p = ProgressBarTrajectory(traj, style='tqdm')
        pt.molsurf(p)

        try:
            p = ProgressBarTrajectory(traj, style=tqdm_notebook)
            pt.molsurf(p)
        except traitlets.traitlets.TraitError:
            pass

    # make sure not loading all coordinates from TrajectoryIterator

    traj2 = pt.iterload(fn('tz2.nc'), fn('tz2.parm7'))
    traj2._size_limit_in_GB = traj2._estimated_GB - 0.001

    with pytest.raises(MemoryError):
        traj2.xyz

    p2 = ProgressBarTrajectory(traj2)
    aa_eq(pt.rmsd(p2), pt.rmsd(traj2))


@unittest.skipUnless(has_ipython, 'require IPython')
def test_pickle_progress_bar():
    traj = pt.datafiles.load_tz2()

    with tempfolder():
        t2 = ProgressBarTrajectory(traj)
        pt.to_pickle(t2, 'test.pk')
        t3 = pt.read_pickle('test.pk')
        data0 = pt.rmsd(traj)
        data1 = pt.rmsd(t2)
        data2 = pt.rmsd(t3)
        aa_eq(data0, data1)
        aa_eq(data0, data2)
