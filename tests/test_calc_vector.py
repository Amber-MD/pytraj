from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir

class Test(unittest.TestCase):
    def test_0(self):
        import pytraj.common_actions as pyca
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        v0 = pyca.calc_vector("@CA @N,C,O", traj)
        print (v0.to_ndarray())
        v1 = traj.calc_vector("@CA @N,C,O")
        print (v1.to_ndarray().shape)
        print (v1.tolist())

if __name__ == "__main__":
    unittest.main()
