#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
# from pytraj.testing import aa_eq
from pytraj.testing import cpptraj_test_dir


class TestTI(unittest.TestCase):
    def test_ti(self):
        dvdl_fn = cpptraj_test_dir + '/Test_TI/dvdl.dat'
        options = 'nq 9'
        pt.ti(dvdl_fn, options)


if __name__ == "__main__":
    unittest.main()
