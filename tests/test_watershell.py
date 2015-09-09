import os
import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj import adict


class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload(".//data/tz2.truncoct.nc",
                             ".//data/tz2.truncoct.parm7")
        dslist = DatasetList()
        dflist = DataFileList()
        adict['watershell'](current_frame=traj,
                            command="!:WAT out ./output/_ws.agr",
                            dslist=dslist,
                            dflist=dflist)

    def test_1(self):
        traj = mdio.iterload(".//data/tz2.truncoct.nc",
                             ".//data/tz2.truncoct.parm7")
        from pytraj.common_actions import calc_watershell

        d0 = calc_watershell(traj, '!:WAT')


if __name__ == "__main__":
    unittest.main()
