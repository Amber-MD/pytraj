#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestTools(unittest.TestCase):

    def test_as_2darray(self):
        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")
        aa_eq(pt.tools.as_2darray(traj), pt.tools.as_2darray(traj.xyz))

        assert pt.tools.as_2darray(traj).ndim == 2, 'ndim must be 2'


if __name__ == "__main__":
    unittest.main()
