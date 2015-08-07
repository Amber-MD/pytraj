import unittest
from pytraj.compat import izip
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):
    def test_0(self):
        # seriously I need to do this test
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        trajsaved = mdio.iterload(
            "./data/fit_to_1stframe.Tc5b.x", "./data/Tc5b.top")

        for frame in trajsaved:
            print(frame[0])

        for f0, f1 in izip(traj, trajsaved):
            print(f0[0], f1[0])


if __name__ == "__main__":
    unittest.main()
