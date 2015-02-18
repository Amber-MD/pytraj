import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        ref = mdio.load("./data/Tc5b.crd", "./data/Tc5b.top")[0]
        traj.top.set_reference_coord(ref)
        act = adict['rmsd']
        dslist = DataSetList()
        act("rms reference :3-18@CA", traj, dslist=dslist)
        print ((dslist[0][:]))

        act = adict['rmsd']
        dslist = DataSetList()
        act("rms first :3-18@CA", traj, dslist=dslist)
        print ((dslist[0][:]))

if __name__ == "__main__":
    unittest.main()
