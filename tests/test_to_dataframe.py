from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.utils import has_
from pytraj.testing import test_if_having


class Test(unittest.TestCase):
    @test_if_having("pandas")
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        act = adict['multidihedral']
        dslist = DatasetList()
        act("phi", traj, dslist=dslist)

        dframe = dslist.to_dataframe()

        self.assertRaises(
            TypeError, lambda: dslist.to_dataframe(engine='xray'))


if __name__ == "__main__":
    unittest.main()
