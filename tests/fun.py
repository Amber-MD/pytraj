#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestFun(unittest.TestCase):
    def test_fun(self):
        '''for anything that does not need serious testing
        '''
        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")

        state = pt.load_pipeline(traj, '''distance :3 :2''').compute()

if __name__ == "__main__":
    unittest.main()
