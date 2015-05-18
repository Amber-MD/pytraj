from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir, duplicate_traj
import pytraj.common_actions as pyca

class Test(unittest.TestCase):
    @test_if_having("numpy")
    def test_0(self):
        from pytraj._xyz import XYZ
        import numpy as np
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        def test_assgin():
            traj.xyz[0] += 1.

        self.assertRaises(NotImplementedError, lambda: test_assgin())
        arr = np.asarray(traj.xyz)
        assert isinstance(arr, np.ndarray)
        xyz = traj.xyz
        assert xyz.shape == (traj.n_frames, traj.n_atoms, 3)
        assert isinstance(xyz, XYZ) == True
        assert isinstance(xyz.tolist(), list) == True

        # try to make array
        arr = np.asarray(xyz)

        # try to do some basic math
        xyz/3.

        def try_to_set_atts():
            traj.xyz = traj.xyz + 1.

        def try_to_setitem():
            traj.xyz[0, 0, 0] = 1.
        self.assertRaises(NotImplementedError, lambda: try_to_setitem())

if __name__ == "__main__":
    unittest.main()
