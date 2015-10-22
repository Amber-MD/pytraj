#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq

class TestIterFrame(unittest.TestCase):
    def test_iterframe(self):
        '''test iterframe for both Trajectory and TrajectoryIterator
        '''

        orig_traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")
        # iterframe (already in doctest), just throwing raise to increase coverage score

        for traj in [orig_traj, orig_traj[:]]:
            self.assertRaises(ValueError, lambda: traj.iterframe(rmsfit='crazy'))


if __name__ == "__main__":
    unittest.main()
    # nosetests --with-coverage --cover-package pytraj -vs .
    # nosetests -vs --processes 6 --process-timeout 200 .
    # nosetests -vs --processes 6 --process-timeout 200 --with-coverage --cover-package pytraj .
