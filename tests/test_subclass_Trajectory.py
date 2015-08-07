from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca
from pytraj._action_in_traj import ActionTrajectory


class Test(unittest.TestCase):
    def test_0(self):
        from pytraj import Trajectory
        # Aim: no segmentation fault
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        class SimpleTrajetory(Trajectory, ActionTrajectory):
            pass

        straj = SimpleTrajetory()
        straj.top = traj.top.copy()
        straj.load(traj)
        aa_eq(straj.xyz, traj.xyz)

        # make sure that when calling `load`, we make a copy
        straj[0, 0, 0] = 100.
        assert straj[0, 0, 0] == 100.
        assert traj[0, 0, 0] != 100.


if __name__ == "__main__":
    unittest.main()
