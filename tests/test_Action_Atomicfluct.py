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
        traj = mdio.iterload("./data/tz2.nc", "./data/tz2.parm7")
        d0 = pyca.calc_atomicfluct(traj, "byres @C,CA,N bfactor")
        print(d0)
        print(d0[0].cpptraj_dtype)
        assert (d0[0].cpptraj_dtype == 'xymesh')
        print(d0.size)
        print(d0.tolist())


if __name__ == "__main__":
    unittest.main()
