import os
import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj import adict

class Test(unittest.TestCase):
    def test_0(self):
        path = "../NOGIT/CpptrajTest"
        if os.path.exists("../NOGIT/CpptrajTest"):
            print (path)
            traj = mdio.load("../NOGIT/CpptrajTest/tz2.truncoct.nc", 
                             "../NOGIT/CpptrajTest/tz2.truncoct.parm7")
            dslist = DataSetList()
            dflist = DataFileList()
            adict['watershell']("!:WAT out ./output/_ws.agr", traj, 
                                dslist=dslist, dflist=dflist)
            print (dslist[0][:])
            print (dslist[1][:])
            print (traj.top[":128"])
            dfile = dflist.get_datafile('./output/_ws.agr')
            print (dir(dfile))
            print (dfile.dtype)
            dflist.write_all_datafiles()

if __name__ == "__main__":
    unittest.main()
