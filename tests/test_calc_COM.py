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
    @test_if_having("numpy")
    def test_0(self):
        import numpy as np

        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        d0 = pyca.calc_COM(traj)
        d1 = pyca.calc_center_of_mass(traj)
        print (d0)
        arr = d0.to_ndarray()
        print (arr)

        for frame in traj:
            print (frame.VGeometricCenter(traj.top("*")).tolist())

        saved_d0 = np.loadtxt("./data/vec.out", skiprows=1, usecols=(1, 2, 3))
        print (saved_d0)

        assert_almost_equal(arr.flatten(), saved_d0.flatten())
        assert_almost_equal(d1.to_ndarray().flatten(), saved_d0.flatten())

if __name__ == "__main__":
    unittest.main()
