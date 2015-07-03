import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj import calculate
from pytraj import *


class Test(unittest.TestCase):

    def test_0(self):
        from pytraj.datasets.DataSetList import DataSetList
        traj = mdio.iterload("./data/tz2.crd", "./data/tz2.parm7")[:]
        # NOT SURE CORRECTLY YET
        # FIXME: Segmentation fault (core dumped)
        #print (d0)
        act = adict['pairdist']
        dslist = DataSetList()
        act('mask "*" mask2 "*" delta 0.1 out ./output/pairdist.dat',
            traj, dslist=dslist)
        # act.print_output()
        print(dslist.size)
        #print (dslist[0])

if __name__ == "__main__":
    unittest.main()
