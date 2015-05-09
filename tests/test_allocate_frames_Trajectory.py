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
        from pytraj import Trajectory
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        traj2 = Trajectory()
        traj2._allocate(traj.n_frames, traj.n_atoms)
        assert (traj2.shape == traj.shape)

        traj2.top = traj.top.copy()
        traj2.update_xyz(traj.xyz)
        aa_eq(traj2.xyz, traj.xyz)

if __name__ == "__main__":
    unittest.main()
