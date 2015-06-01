from __future__ import print_function
import unittest # pragma no_test for travis
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir, duplicate_traj
import pytraj.common_actions as pyca

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        xyz = traj.xyz[:]
        # x-coords for all atoms
        # should it?
        #aa_eq(traj[:, :, 0], xyz[:, :, 0])

if __name__ == "__main__":
    unittest.main()
