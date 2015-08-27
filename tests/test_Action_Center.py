from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir, duplicate_traj
import pytraj.common_actions as pyca
from pytraj.compat import zip


class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        fa = traj[:]
        # load 2 frames
        saved_traj = mdio.load("./data/tz2.center_mass.nc", traj.top)
        fa.center(":1 mass")

        aa_eq(fa[:2].xyz, saved_traj.xyz, decimal=5)


if __name__ == "__main__":
    unittest.main()
