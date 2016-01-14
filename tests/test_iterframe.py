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
            self.assertRaises(ValueError,
                              lambda: traj.iterframe(rmsfit='crazy'))

            # rmsfit is an int
            t0 = orig_traj[:].rmsfit(ref=3)
            aa_eq(
                pt.rmsd_nofit(
                    traj(rmsfit=3),
                    ref=orig_traj[-1]),
                pt.rmsd_nofit(t0, ref=orig_traj[-1]))

        # test TypeError if not has n_frames info
        t0 = orig_traj[:]

        def int_gen(k):
            for i in range(k):
                yield i

        fi = pt.iterframe(t0, frame_indices=int_gen(3))
        aa_eq(pt.radgyr(fi, top=traj.top), pt.radgyr(orig_traj[:3]))


if __name__ == "__main__":
    unittest.main()
