import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        farray = traj[:]
        print (traj[0, 0])
        print (traj[-1, 0])

        print (traj)
        for frame in traj.frame_iter():
            print (frame[0])

        for frame0 in farray.frame_iter():
            print (frame0[0])

if __name__ == "__main__":
    unittest.main()
