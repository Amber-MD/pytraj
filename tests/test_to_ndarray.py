from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj.exceptions import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal, eq
from pytraj.decorators import no_test, test_if_having
from pytraj.utils import _import_numpy

has_np, np = _import_numpy()

class Test(unittest.TestCase):
    @test_if_having("numpy")
    def test_0(self):
        import numpy as np
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        # frame
        frame0 = traj[0]
        arr0 = frame0.to_ndarray()
        assert isinstance(arr0, np.ndarray) == True
        assert_almost_equal(arr0.flatten(), frame0.coords)

        # test memview
        arr0[0, 0] = 1000.
        assert arr0[0, 0] == frame0[0, 0] == 1000.

    def test_1(self):
        from pytraj.common_actions import calc_multidihedral
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        d0 = calc_multidihedral("psi", traj)
        print (d0)
        if not has_np:
            self.assertRaises(PytrajConvertError, lambda: dslist.to_ndarray())


if __name__ == "__main__":
    unittest.main()
