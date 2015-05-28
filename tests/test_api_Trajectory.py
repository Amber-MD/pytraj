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
        api_traj = api.Trajectory(traj)

        # test xyz
        aa_eq(api_traj.xyz, traj.xyz)

        # test object lifetime
        aa_eq(api_traj[0].xyz, api_traj.xyz[0])

        # test Box
        assert (api_traj.has_box() == True)
        boxes = traj.box_to_ndarray()
        for i, frame in enumerate(api_traj):
            assert (frame.has_box() == True)
            f_blist = frame.box.tolist()
            aa_eq(f_blist, boxes[i].tolist())

        # test join
        api_traj += traj
        assert api_traj.n_frames == traj.n_frames * 2
        aa_eq(api_traj.xyz[traj.n_frames:], api_traj.xyz[:traj.n_frames])

        # test autoimage
        # make Trajectory from TrajectoryIterator
        fa = traj[:]
        t_api2 = api.Trajectory(fa)
        t_api2.autoimage()
        fa.autoimage()
        saved_traj = mdio.iterload("./data/tz2.autoimage.nc", "./data/tz2.ortho.parm7")
        for f1, f2 in zip(fa, t_api2):
            aa_eq(f1.box.tolist(), f2.box.tolist())
            print (f1.rmsd(f2))
        aa_eq(t_api2.xyz, fa.xyz, decimal=1)

        # make sure to reproduce cpptraj's output too
        aa_eq(saved_traj.xyz, fa.xyz)

        # test slicing to get correct Box
        t_api2_sliced = t_api2[:10]
        aa_eq(t_api2_sliced._boxes, t_api2._boxes[:10])

if __name__ == "__main__":
    unittest.main()
