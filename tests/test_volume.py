#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestVolume(unittest.TestCase):

    def test_volume(self):
        traj = pt.iterload("data/tz2.ortho.nc", "data/tz2.ortho.parm7")
        state = pt.load_cpptraj_state('''
        volume''', traj)
        state.run()

        vol = pt.volume(traj)
        aa_eq(state.data[-1], vol)


if __name__ == "__main__":
    unittest.main()
