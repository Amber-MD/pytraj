from __future__ import print_function
import unittest
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.testing import test_if_having
from pytraj import io as mdio


class Test(unittest.TestCase):
    @test_if_having("numpy")
    def test_0(self):
        import numpy as np
        traj = mdio.iterload("./data/tz2.nc", "./data/tz2.parm7")
        top = traj.top

        # test for last frame
        top.set_reference_frame(traj[-1])
        indices = top.select(":3@CA <:3.0").indices

        saved_indices = np.loadtxt(
            "./data/mask.tz2.dat",
            skiprows=1,
            usecols=(1, ))
        # subtract by '1' since cpptraj uses "1" as starting index for output
        saved_indices = saved_indices - 1
        assert_almost_equal(indices, saved_indices)


if __name__ == "__main__":
    unittest.main()
