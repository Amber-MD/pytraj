from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        print ("test_0")
        traj = mdio.load("./data/tz2.truncoct.nc", "./data/tz2.truncoct.parm7")
        print (traj.size)
        act = adict['vector']
        dslist = DataSetList()
        act("v0 principal x @CA", traj, dslist=dslist)
        print ('dslist.size = ', dslist.size)
        print (dslist.get_dtypes())

        for d0 in dslist:
            for vec in d0.data:
                print (vec.tolist())

    def test_1(self):
        print ("test_1")
        traj = mdio.load("./data/tz2.truncoct.nc", "./data/tz2.truncoct.parm7")
        print (traj.size)
        act = adict['vector']
        dslist = DataSetList()
        act("v0 principal x @CA", traj, dslist=dslist)
        act("v1 principal y @CA", traj, dslist=dslist)
        act("v2 principal z @CA", traj, dslist=dslist)
        print ('dslist.size = ', dslist.size)

        for d0 in dslist:
            for vec in d0.data:
                print (vec.tolist())

        print (dslist['v2'][0][:])

if __name__ == "__main__":
    unittest.main()
