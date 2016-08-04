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

        command = "doorder refdens 0.033422885325 gridcntr 1.5 1.0 0.0 griddim 34 44 36 gridspacn 0.50"
        pt.all_actions.gist(traj, command)

        def get_gist_data(filename):
            try:
                return np.loadtxt(filename, skiprows=7).transpose()
            except ValueError:
                print(filename)

        gist_outputs = [
                'gist-gH.dx',
                'gist-gO.dx',
                'gist-neighbor-norm.dx',
                'gist-order-norm.dx',
        ]

        for filename in gist_outputs:
            saved_gist_gH_filename = cpptraj_test_dir + '/Test_GIST/' + filename + '.save'
            saved_data = get_gist_data(saved_gist_gH_filename)

            output = get_gist_data(filename)

            aa_eq(output, saved_data)

if __name__ == "__main__":
    unittest.main()
