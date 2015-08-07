from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.testing import test_if_having
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):
    @test_if_having("numpy")
    def test_0(self):
        import numpy as np
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        frame = traj[0]
        assert_almost_equal(frame.coords, np.array(frame.tolist()).flatten())
        assert len(frame.tolist()) == frame.n_atoms


if __name__ == "__main__":
    unittest.main()
