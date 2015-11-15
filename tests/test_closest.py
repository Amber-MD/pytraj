#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestClosest(unittest.TestCase):
    def test_closest(self):
        # raise if not has solvent
        traj0 = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")
        self.assertRaises(RuntimeError, lambda: pt.closest(traj0))

        traj = pt.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        fi, top = pt.closest(traj)
        for frame in fi:
            assert isinstance(frame, pt.Frame), 'must be Frame'

if __name__ == "__main__":
    unittest.main()
