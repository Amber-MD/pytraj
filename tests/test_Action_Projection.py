from __future__ import print_function

# turn off test: don't know how to test yet
#import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        act = adict['projection']
        dslist = DataSetList()
        act("@CA", traj, dslist=dslist)
        print(dslist.size)
        print(dslist.get_legends())


if __name__ == "__main__":
    unittest.main()
