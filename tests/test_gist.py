#!/usr/bin/env python

from __future__ import print_function
import unittest
import numpy as np
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.testing import cpptraj_test_dir


class TestGist(unittest.TestCase):

    def test_gist(self):
        traj = pt.iterload("data/tz2.ortho.nc", "data/tz2.ortho.parm7", frame_slice=(0, 10))
        traj.autoimage('origin')

        command = "doorder doeij refdens 0.033422885325 gridcntr 1.44 0.67 0.29 \
                     griddim 10 12 10 gridspacn 2.0"
        data = pt.all_actions.gist(traj, command)

if __name__ == "__main__":
    unittest.main()
