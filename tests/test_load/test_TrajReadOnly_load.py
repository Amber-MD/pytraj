import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):

    def test_0(self):
        traj = TrajectoryIterator("./data/md1_prod.Tc5b.x",
                                  "./data/Tc5b.top",
                                  frame_slice=(1, 8, 2))


if __name__ == "__main__":
    unittest.main()
