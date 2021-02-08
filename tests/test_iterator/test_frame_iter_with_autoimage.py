from __future__ import print_function
import unittest
from pytraj import *
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq


class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload(fn('tz2.ortho.nc'), fn('tz2.ortho.parm7'))

        fa1 = traj[:]

        # do inplace-autoimage for Trajectory
        fa1.autoimage()

        # fa2 will store frame when iterating `traj` (TrajectoryIterator)
        fa2 = Trajectory()
        fa2.top = traj.top.copy()

        # frame_iter
        for frame in traj(autoimage=True):
            fa2.append(frame)

        aa_eq(fa1.xyz, fa2.xyz)

        # frame_iter with mask and autoimage
        fa3 = fa2.copy()
        # fa3 is a copy of autoimaged fa2. then we strip all but CA atoms
        # just want to make sure we can use `mask`
        fa3.strip("!@CA")
        fa4 = Trajectory()
        fa4.top = fa3.top.copy()

        for frame in traj(mask='@CA', autoimage=True):
            fa4.append(frame)

        aa_eq(fa3.xyz, fa4.xyz)

        # frame_iter with mask and autoimage, and rmsfit
        # fa3 is a copy of autoimaged fa2. then we strip all but CA atoms
        # just want to make sure we can use `mask`

        fa3 = traj[:]
        fa3.autoimage()
        fa3.rmsfit(ref=5, mask='@CB')
        fa3.strip("!@CA")

        fa4 = Trajectory()
        fa4.top = fa3.top.copy()

        ref1 = traj[5].copy()
        for frame in traj(mask='@CA', autoimage=True, rmsfit=(ref1, '@CB')):
            fa4.append(frame)

        aa_eq(fa3.xyz, fa4.xyz)
        for f0, f1 in zip(fa3, fa4):
            assert f0.rmsd_nofit(f1) < 1E-7


if __name__ == "__main__":
    unittest.main()
