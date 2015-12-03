#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class Test(unittest.TestCase):

    def test_0(self):
        traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")


if __name__ == "__main__":
    unittest.main()
