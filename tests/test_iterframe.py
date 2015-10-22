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

            # rmsfit is an int
            t0 = orig_traj[:].rmsfit(3)
            aa_eq(pt.rmsd_nofit(traj(rmsfit=3), orig_traj[-1]),
                  pt.rmsd_nofit(t0, orig_traj[-1]))



if __name__ == "__main__":
    unittest.main()
    # nosetests --with-coverage --cover-package pytraj -vs .
    # nosetests -vs --processes 6 --process-timeout 200 .
    # nosetests -vs --processes 6 --process-timeout 200 --with-coverage --cover-package pytraj .
