from __future__ import print_function
import pytraj as pt
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca
from pytraj.compat import zip


class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        print(traj)

        fa1 = traj[:]

        # do inplace-autoimage for Trajectory
        fa1.autoimage()

        # fa2 will store frame when iterating `traj` (TrajectoryIterator)
        fa2 = Trajectory()
        fa2.top = traj.top.copy()

        # frame_iter
        for frame in traj(autoimage=True):
            fa2.append(frame, copy=True)

        aa_eq(fa1.xyz, fa2.xyz)

        # frame_iter with mask and autoimage
        fa3 = fa2.copy()
        # fa3 is a copy of autoimaged fa2. then we strip all but CA atoms
        # just want to make sure we can use `mask`
        fa3.strip_atoms("!@CA")
        fa4 = Trajectory()
        fa4.top = fa3.top.copy()

        for frame in traj(mask='@CA', autoimage=True):
            fa4.append(frame, copy=True)

        print(fa3, fa4)
        print(fa3.top, fa4.top)
        print(fa3.shape, fa4.shape)
        aa_eq(fa3.xyz, fa4.xyz)

        # frame_iter with mask and autoimage, and rmsfit
        # fa3 is a copy of autoimaged fa2. then we strip all but CA atoms
        # just want to make sure we can use `mask`

        fa3 = traj[:]
        fa3.autoimage()
        fa3.rmsfit(5, '@CB')
        fa3.strip_atoms("!@CA")

        fa4 = Trajectory()
        fa4.top = fa3.top.copy()

        ref1 = traj[5].copy()
        for frame in traj(mask='@CA', autoimage=True, rmsfit=(ref1, '@CB')):
            fa4.append(frame, copy=True)

        print(fa3, fa4)
        aa_eq(fa3.xyz, fa4.xyz)
        for f0, f1 in zip(fa3, fa4):
            assert f0.rmsd_nofit(f1) < 1E-7
            print(f0.rmsd_nofit(f1))


if __name__ == "__main__":
    unittest.main()
