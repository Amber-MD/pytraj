#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.testing import cpptraj_test_dir


class TestForce(unittest.TestCase):

    def test_nosegfault_for_force(self):
        fn = cpptraj_test_dir + '/Test_systemVF/systemVF.nc'
        tn = cpptraj_test_dir + '/Test_systemVF/systemVF.parm7'
        traj = pt.iterload(fn, tn)

        for f in traj:
            assert f.has_force(), 'must have force'


if __name__ == "__main__":
    unittest.main()
