from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import has_
from pytraj.utils.check_and_assert import assert_almost_equal, eq
from pytraj.testing import aa_eq


class Test(unittest.TestCase):

    def test_0(self):
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        if has_("numpy"):
            import numpy as np
            from numpy.testing import assert_almost_equal
            arr0 = traj.xyz
            aa_eq(arr0, traj[:, :, :].xyz)

            # create Trajectory
            farray = traj[:]
            aa_eq(arr0, farray[:, :, :].xyz)
        else:
            pass


if __name__ == "__main__":
    unittest.main()
