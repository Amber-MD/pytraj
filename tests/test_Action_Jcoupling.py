import os
import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj import adict
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj import DataSetList, DataFileList

class Test(unittest.TestCase):
    def test_0(self):
        try:
            os.environ['AMBERHOME']
            has_amberhome = True
        except:
            has_amberhome = False
        if has_amberhome:
            print ("has_amberhome, doing calculation")
            traj = mdio.load("./data/tz2.nc", "./data/tz2.parm7")
            frame = traj[0]
            dslist = DataSetList()
            dflist = DataFileList()
            d0 = adict['jcoupling']("out ./output/test_jcoupling.out",
                                   frame, traj.top, 
                                   dslist=dslist, dflist=dflist)
            print (dslist.size)
            print (dslist[0].data)
            print (dslist[80])
            print (dslist.get_legends())
            print (dslist['LYS:8_C-CA-CB-HB2'][0][:])

            # another way
            dslist = DataSetList()
            from pytraj.common_actions import calc_jcoupling
            d0 = calc_jcoupling("out ./output/test_jcoupling_1.out", traj, dslist=dslist)
            print (d0.legend)
            print (dslist.size)


if __name__ == "__main__":
    unittest.main()
