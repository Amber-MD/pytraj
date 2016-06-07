#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestStrippedTrajectory(unittest.TestCase):

    def test_stripped_trajectory(self):
        traj = pt.iterload("data/tz2.ortho.nc", "data/tz2.ortho.parm7")
        straj = traj.autoimage().superpose('@CA').strip(":WAT")

        aa_eq(traj['!:WAT'].xyz, straj.xyz)

if __name__ == "__main__":
    unittest.main()
