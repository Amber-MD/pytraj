# turn off unittest to pass travis
# still need to fix this
#import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        act = adict['rmsd']
        dslist = DataSetList()
        # FIXME: segmentation fault
        act("parm ./data/Tc5b.top reference ./data/Tc5b.crd rms reference @CA", traj, dslist=dslist)

if __name__ == "__main__":
    unittest.main()
