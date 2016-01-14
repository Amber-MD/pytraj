#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestPipeline(unittest.TestCase):

    def test_pieline(self):
        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")

        fi = pt.pipe(traj, ['autoimage', ])
        aa_eq(pt.get_coordinates(fi), traj[:].autoimage().xyz)

        fi = pt.pipe(traj, ['autoimage', ], frame_indices=[3, 5])
        aa_eq(pt.get_coordinates(fi), traj[[3, 5]].autoimage().xyz)


if __name__ == "__main__":
    unittest.main()
