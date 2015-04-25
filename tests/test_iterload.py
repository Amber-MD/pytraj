from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca
from pytraj.six_2 import izip

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        itertraj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top", 0, -1, 1, "")

        for idx, (f0, f1) in enumerate(izip(traj, itertraj)):
            assert_almost_equal(f0.coords, f1.coords)
        assert idx == traj.n_frames - 1

if __name__ == "__main__":
    unittest.main()
