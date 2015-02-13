import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")[:]

        frame = traj[0]
        print (frame[0, :])

        for arr0 in frame:
            arr0[0] = 0.0

        # make sure that we update frame too
        for i in range(frame.n_atoms):
            assert (frame[i, 0]) == 0.0

        print ("test frame_iter")
        for f0 in frame.frame_iter():
            print (f0)

if __name__ == "__main__":
    unittest.main()
