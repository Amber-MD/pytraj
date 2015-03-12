import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        # TODO: add assertion
        # FIXME: results seem wrong
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        act = adict['dihedralscan']
        dslist = DataSetList()
        act("@CA out ./output/_test_Dihscan.dat", traj, dslist=dslist)
        print (dslist.size)
        d0 = dslist[0]
        print (d0)
        print (d0[:])
        act.print_output()

if __name__ == "__main__":
    unittest.main()
