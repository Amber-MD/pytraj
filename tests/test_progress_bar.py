#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils.progress import ProgressBarTrajectory
from pytraj.utils import aa_eq
import pytest

try:
    import tqdm
    from tqdm import tqdm_notebook
except ImportError:
    tqdm = None


class TestProgressLog(unittest.TestCase):

    def test_progress_log(self):
        """test_progress_log: simple test, just to make sure it is runnable
        """
        traj = pt.datafiles.load_tz2()

        p = ProgressBarTrajectory(traj, style='basic', every=20)
        pt.molsurf(p)

        p = ProgressBarTrajectory(traj, style='bar', every=20)
        pt.molsurf(p)

        p = ProgressBarTrajectory(traj, style='circle', every=20)
        pt.molsurf(p)

        # need to run on Jupyter notebook
        if tqdm is not None:
            p = ProgressBarTrajectory(traj, style='tqdm')
            pt.molsurf(p)

            p = ProgressBarTrajectory(traj, style=tqdm_notebook)
            pt.molsurf(p)

        # make sure not loading all coordinates from TrajectoryIterator

        traj2 = pt.iterload('data/tz2.nc', 'data/tz2.parm7')
        traj2._size_limit_in_GB = traj2._estimated_GB - 0.001

        
        with pytest.raises(MemoryError):
            traj2.xyz

        p2 = ProgressBarTrajectory(traj2)
        aa_eq(pt.rmsd(p2), pt.rmsd(traj2))
        


if __name__ == "__main__":
    unittest.main()
