import unittest
from pytraj import adict
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")[:]
        traj2 = traj.copy()
        print (traj['@CA'][0, 0])

        # test action without expliti specifying Topology
        adict['rotate']('@CA x 60 y 60 z 60', traj)
        print (traj['@CA'][0, 0])
        adict['rotate']('@CA x 60 y 60 z 60', traj2, traj2.top)
        assert_almost_equal(traj2['@CA'][0, 0], traj['@CA'][0, 0])

if __name__ == "__main__":
    unittest.main()
