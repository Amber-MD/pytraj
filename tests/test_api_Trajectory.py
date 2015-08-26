from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca
from pytraj import api
from pytraj.compat import zip


class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        api_traj = traj[:]

        # test xyz
        aa_eq(api_traj.xyz, traj.xyz)

        # test object lifetime
        aa_eq(api_traj[0].xyz, api_traj.xyz[0])

        # test Box
        assert (api_traj.has_box() == True)
        boxes = traj.unitcells
        for i, frame in enumerate(api_traj):
            assert (frame.has_box() == True)
            f_blist = frame.box.tolist()
            aa_eq(f_blist, boxes[i].tolist())

        # test autoimage
        # make Trajectory from TrajectoryIterator
        fa = traj[:]
        fa.autoimage()
        saved_traj = mdio.iterload(
            "./data/tz2.autoimage.nc", "./data/tz2.ortho.parm7")

        # make sure to reproduce cpptraj's output too
        aa_eq(saved_traj.xyz, fa.xyz)


if __name__ == "__main__":
    unittest.main()
