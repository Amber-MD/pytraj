import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        # TODO: add assert
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist = DataSetList()
        act = adict['pairwise']
        act("@CA", traj, dslist=dslist)
        print (dslist.size)

        for ds in dslist:
            print (ds)

        d3 = dslist[3]
        print (d3.dtype)
        #$for i in range(d3.size):
        #$    print (d3[i])
        #$#print (dslist[3].data)

if __name__ == "__main__":
    unittest.main()
