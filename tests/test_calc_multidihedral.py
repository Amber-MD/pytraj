from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):

    def test_0(self):
        from pytraj import DatasetList
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        from pytraj.all_actions import calc_multidihedral
        command = "resrange 2-19 phi psi"
        d0 = calc_multidihedral(traj, command)
        d1 = calc_multidihedral(traj, command)
        assert isinstance(d0, DatasetList) == True
        assert (len(d0.keys()) == len(d1.keys()))

        d3 = calc_multidihedral(traj)

        import numpy as np
        d0np = calc_multidihedral(traj, command, dtype='ndarray')
        self.assertIsInstance(d0np, np.ndarray)


if __name__ == "__main__":
    unittest.main()
