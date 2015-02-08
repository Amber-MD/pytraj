import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj import adict
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj import DataSetList, DataFileList

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/tz2.nc", "./data/tz2.parm7")
        frame = traj[0]
        dslist = DataSetList()
        dflist = DataFileList()
        d0 = adict['jcoupling']("out ./output/test_jcoupling.out", 
                               frame, traj.top, 
                               dslist=dslist, dflist=dflist)
        print (dslist[0])
        dflist.write_all_datafiles()

if __name__ == "__main__":
    unittest.main()
