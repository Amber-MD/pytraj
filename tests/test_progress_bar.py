#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils.progress import ProgressBarTrajectory
from pytraj.utils import eq, aa_eq
from tqdm import tqdm_notebook


class TestProgressLog(unittest.TestCase):

    def test_progress_log(self):
        """test_progress_log: simple test, just to make sure it is runnable
        """
        traj = pt.datafiles.load_tz2()
        p = ProgressBarTrajectory(traj, style='bar', every=20)
        pt.molsurf(p)

        p = ProgressBarTrajectory(traj, style='circle', every=20)
        pt.molsurf(p)

        # need to run on Jupyter notebook
        # p = ProgressBarTrajectory(traj, style='tqdm')
        # pt.molsurf(p)

        # p = ProgressBarTrajectory(traj, style=tqdm_notebook)
        # pt.molsurf(p)


if __name__ == "__main__":
    unittest.main()
