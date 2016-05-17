#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class Test(unittest.TestCase):

    def test_0(self):
        traj = pt.load("./data/DPDP.nc", "data/DPDP.parm7")
        traj.superpose('@C,CA,N', ref=0)
        data = pt.all_actions._wavelet(traj, 'nb 10 s0 2 ds 0.25 type morlet correction 1.01 chival 0.25 :1-22 name DPDP cluster minpoints 66 epsilon 10.0')
        print(data)
        
if __name__ == "__main__":
    unittest.main()
