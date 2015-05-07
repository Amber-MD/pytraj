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
        traj = mdio.load("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        traj.autoimage()
        ds = pyca.search_hbonds(traj)
        d0 = ds[0]
        assert d0.count(4) == 3
        d1 = ds[1]
        assert d1.count(1) == 7

if __name__ == "__main__":
    unittest.main()
