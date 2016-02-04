from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io
from pytraj.utils.check_and_assert import assert_almost_equal

import itertools as it
"""
try not to get segmentation fault error (due to whatever freaking reason)
"""
traj = io.iterload("./data/Tc5b.x", "./data/Tc5b.top")


class Test(unittest.TestCase):

    def test_0_trajiter(self):
        traj = io.load_sample_data("tz2")
        from pytraj.compat import zip

        for idx, (f0, f1) in enumerate(zip(traj, traj)):
            pass
        assert idx > 0

        # tee
        ilist = it.tee(traj, 2)
        for idx, (f0, f1) in enumerate(zip(ilist[0], ilist[1])):
            assert f0.rmsd(f1) < 1E-5

        assert idx == traj.n_frames - 1

        # product
        it.product(traj, repeat=2)


if __name__ == "__main__":
    unittest.main()
