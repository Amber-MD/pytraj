import os
import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj import adict

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load(".//data/tz2.truncoct.nc", 
                         ".//data/tz2.truncoct.parm7")
        dslist = DataSetList()
        dflist = DataFileList()
        adict['watershell'](current_frame=traj, command="!:WAT out ./output/_ws.agr",
                            dslist=dslist, dflist=dflist)
        print (dslist[0][:])
        print (dslist[1][:])

    def test_1(self):
        traj = mdio.load(".//data/tz2.truncoct.nc", 
                         ".//data/tz2.truncoct.parm7")
        from pytraj.common_actions import calc_watershell

        d0 = calc_watershell(traj, '!:WAT')
        print (d0[0][:])
        print (d0[1][:])

if __name__ == "__main__":
    unittest.main()
