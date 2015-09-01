import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class TestTrajinList(TrajinList):
    def __getitem__(self, idx):
        s = 0
        for traj in self:
            if s == idx:
                traj.top = self.top
                return traj
            s += 1


class Test(unittest.TestCase):
    def test_0(self):
        cppstate = CpptrajState()
        cppstate.toplist.add_parm(
            "./data/Test_RemdTraj/ala2.99sb.mbondi2.parm7")
        cppstate.add_trajin("./data/Test_RemdTraj/rem.nc.000")

    def test_1(self):
        trajinlist = TrajinList()
        top = Topology("./data/Test_RemdTraj/ala2.99sb.mbondi2.parm7")
        trajinlist.add_traj("./data/Test_RemdTraj/rem.nc.000", top, "1 4")
        traj = trajinlist[0]
        traj.top = top.copy()

    def test_2(self):
        trajinlist = TestTrajinList()
        top = Topology("./data/Test_RemdTraj/ala2.99sb.mbondi2.parm7")
        trajinlist.add_traj("./data/Test_RemdTraj/rem.nc.000", top, "1 4")
        traj = trajinlist[0]
        for frame in traj:
            pass


if __name__ == "__main__":
    unittest.main()
