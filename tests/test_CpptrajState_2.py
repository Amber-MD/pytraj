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
        #print(cppstate.toplist)
        #print(cppstate.toplist[0])
        # cppstate.add_trajin(ArgList("./data/Test_RemdTraj/rem.nc.000"),
        #                    is_ensemble=True)
        cppstate.add_trajin("./data/Test_RemdTraj/rem.nc.000")
        #print(dir(cppstate))

    def test_1(self):
        trajinlist = TrajinList()
        top = Topology("./data/Test_RemdTraj/ala2.99sb.mbondi2.parm7")
        trajinlist.add_traj("./data/Test_RemdTraj/rem.nc.000", top, "1 4")
        #print('topology for trajinlist = ', trajinlist.top)
        traj = trajinlist[0]
        traj.top = top.copy()
        #print(traj)
        #print(traj[0])

    def test_2(self):
        trajinlist = TestTrajinList()
        #print(dir(trajinlist))
        top = Topology("./data/Test_RemdTraj/ala2.99sb.mbondi2.parm7")
        trajinlist.add_traj("./data/Test_RemdTraj/rem.nc.000", top, "1 4")
        #print(trajinlist.size)
        #print('topology for trajinlist = ', trajinlist.top)
        traj = trajinlist[0]
        #print(traj)
        #print(traj[0])
        for frame in traj:
            #print(frame)


if __name__ == "__main__":
    unittest.main()
