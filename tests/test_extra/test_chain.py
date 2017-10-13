#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq


class TestChainofCommands(unittest.TestCase):
    def test_center_autoimagesuperpose(self):
        traj = pt.iterload(fn('tz2.ortho.nc'), fn('tz2.ortho.parm7'))
        t0 = traj[:]
        t1 = traj[:]

        pt.center(t0)
        pt.autoimage(t0)
        pt.superpose(t0)

        aa_eq(pt.center(t1).autoimage().superpose().xyz, t0.xyz)


if __name__ == "__main__":
    unittest.main()
