from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")

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

        aa_eq(fa3.xyz, fa4.xyz)

if __name__ == "__main__":
    unittest.main()
