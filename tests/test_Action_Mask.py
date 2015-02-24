import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        act = adict['mask']
        dslist = DataSetList()
        dflist = DataFileList()
        act('mask "(:5 <:3.0)" maskout ./output/test_Action_Mask.dat', traj, dslist=dslist, 
           dflist=dflist)
        print (dslist.size)
        dflist.write_all_datafiles()

if __name__ == "__main__":
    unittest.main()
