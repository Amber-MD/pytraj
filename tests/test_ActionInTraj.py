from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca

class Test(unittest.TestCase):
    def test_0(self):
        from pytraj.misc import get_atts
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        print (traj)
        print (get_atts(traj))
        print (traj.search_hbonds().to_dict())

if __name__ == "__main__":
    unittest.main()
