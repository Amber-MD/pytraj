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
        traj2 = traj.to_mutable_trajectory()
        d0 = pyca.calc_COM(traj, dtype='dataset')
        d1 = pyca.calc_center_of_mass(traj, dtype='dataset')
        d2 = pyca.calc_center_of_mass(traj2, dtype='dataset')
        print(d0)
        arr = d0.to_ndarray()
        print(arr)

        for frame in traj:
            print(frame.center_of_geometry(traj.top("*")).tolist())

        saved_d0 = np.loadtxt("./data/vec.out", skiprows=1, usecols=(1, 2, 3))
        print(saved_d0)

        assert_almost_equal(arr.flatten(), saved_d0.flatten())
        assert_almost_equal(d1.to_ndarray().flatten(), saved_d0.flatten())
        assert_almost_equal(d2.to_ndarray().flatten(), saved_d0.flatten())

        assert_almost_equal(
            traj.calc_COG(dtype='ndarray'), traj2.calc_COG(dtype='ndarray'))


if __name__ == "__main__":
    unittest.main()
