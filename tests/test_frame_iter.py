import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        for i, frame in enumerate(traj.frame_iter(1, 5, 2)):
            pass

        assert i == 2

        for i, frame in enumerate(traj.frame_iter(1, 5, 1)):
            pass

        assert i == 4

        for i, frame in enumerate(traj.frame_iter(stop=8)):
            pass

        assert i == 8

        for i, frame in enumerate(traj.frame_iter(start=7, stop=8)):
            print (frame)

        assert i == 1

if __name__ == "__main__":
    unittest.main()
