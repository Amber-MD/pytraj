from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.testing import cpptraj_test_dir

from pytraj.compat import izip


class Test(unittest.TestCase):

    def test_0(self):
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        itertraj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")

        for idx, (f0, f1) in enumerate(izip(traj, itertraj)):
            assert_almost_equal(f0.xyz, f1.xyz)
        assert idx == traj.n_frames - 1


if __name__ == "__main__":
    unittest.main()
