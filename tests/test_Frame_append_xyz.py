from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        f0 = Frame()
        f0.append_xyz(traj[0].buffer2d)
        print(f0.size)


if __name__ == "__main__":
    unittest.main()
