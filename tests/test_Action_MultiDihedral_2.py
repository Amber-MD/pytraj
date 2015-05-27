from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        act = adict['multidihedral']
        dslist = DataSetList()
        dflist = DataFileList()
        act("phi psi resrange 6-9 out ./output/test_muldih.dat", traj[:2], 
            dslist=dslist, dflist=dflist)
        act.print_output()
        print (dslist.get_legends())
        print (dslist['phi:6'])
        d0 = dslist['phi:6']
        print (d0[:])

if __name__ == "__main__":
    unittest.main()
