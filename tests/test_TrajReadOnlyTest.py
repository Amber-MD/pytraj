import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.TrajReadOnlyTest import TrajReadOnlyTest

class Test(unittest.TestCase):
    def test_0(self):
        traj = TrajReadOnlyTest()
        traj.top = Topology("./data/Tc5b.top")
        print (traj.top)
        traj.load("./data/md1_prod.Tc5b.x")
        traj.load("./data/md1_prod.Tc5b.x")
        traj.load("./data/md1_prod.Tc5b.x")
        traj.load("./data/md1_prod.Tc5b.x")
        print (traj.copy())
        traj.join(traj.copy())
        #print (traj.n_frames)
        #print (traj.n_frames)
        for frame in traj:
            print (frame)

        print (traj.n_frames)
        trajcp = traj.copy()
        print (trajcp.n_frames)

if __name__ == "__main__":
    unittest.main()
