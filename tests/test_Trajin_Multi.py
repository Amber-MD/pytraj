# skip this test since we don't know how to use it
#import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils import assert_almost_equal
from pytraj.trajs.Trajin_Multi import Trajin_Multi

class Test(unittest.TestCase):
    def test_0(self):
        traj = Trajin_Multi()
        traj.top = Topology("./data/Test_RemdTraj/ala2.99sb.mbondi2.parm7")
        traj.load("./data/Test_RemdTraj/rem.nc.000", traj.top, ArgList("ensemble"))
        print (traj)
        traj.ensemble_setup()
        print (traj.ensemble_position(0))
        print (traj.ensemble_position(1))
        print (traj.final_crd_indices())
        print (traj.final_crd_indices())
        print (traj.final_crd_indices())
        print (traj.target_mode)
        farray = FrameArray()
        # FIXME: pytraj has wrong implementation. Check cpptraj code
        traj.get_next_ensemble(farray)
        print (farray)

    def test_1(self):
        traj = TrajinList()
        top = Topology("./data/Test_RemdTraj/ala2.99sb.mbondi2.parm7")
        traj.add_ensemble("./data/Test_RemdTraj/rem.nc.000", top, ArgList("ensemble"))
        print (traj)
        print (traj.mode)
        print (traj.size)
        print (traj.top)
        # FIXME: Segmentation fault (core dumped)
        #print (traj[0])

if __name__ == "__main__":
    unittest.main()
