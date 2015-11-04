from __future__ import print_function
import pytraj as pt
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        import numpy as np

        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        traj2 = traj[:]
        d0 = pyca.calc_COM(traj, dtype='dataset')
        d1 = pyca.calc_center_of_mass(traj, dtype='dataset')
        d2 = pyca.calc_center_of_mass(traj2, dtype='dataset')
        arr = d0.to_ndarray()

        for frame in traj:
            pass

        saved_d0 = np.loadtxt("./data/vec.out", skiprows=1, usecols=(1, 2, 3))

        assert_almost_equal(arr.flatten(), saved_d0.flatten())
        assert_almost_equal(d1.to_ndarray().flatten(), saved_d0.flatten())
        assert_almost_equal(d2.to_ndarray().flatten(), saved_d0.flatten())

        assert_almost_equal(
            pt.center_of_geometry(traj,
                                  dtype='ndarray'),
            pt.center_of_geometry(traj2,
                                  dtype='ndarray'))


if __name__ == "__main__":
    unittest.main()
