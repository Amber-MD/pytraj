import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj import calculate
from pytraj import *

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/tz2.crd", "./data/tz2.parm7")[:]
        # NOT SURE CORRECTLY YET
        # FIXME: Segmentation fault (core dumped)
        #d0 = calculate('pairdist', 'mask "*" delta 0.1', traj)
        #print (d0)
        #act = adict['pairdist']
        #dslist = DataSetList()
        #act.run('mask "*" mask1 "*" delta 0.1', traj, dslist=dslist)
        #print (dslist.size)
        #print (dslist[0])

if __name__ == "__main__":
    unittest.main()
