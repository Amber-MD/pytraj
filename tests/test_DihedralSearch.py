import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        act = adict['dihedralscan']
        dslist = DataSetList()
        act("phi :2-19 psi :2-19", traj, dslist=dslist)
        #print(act.n_frames)
        #print(dslist.size)
        #print(dslist[0])
        #print(dslist[0][:])


if __name__ == "__main__":
    unittest.main()
