#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
import numpy as np


class TestPairDist(unittest.TestCase):
    def test_general(self):
        traj = pt.iterload("./data/tz2.crd", "./data/tz2.parm7")

        data = pt.common_actions.calc_pairdist(traj, delta=0.1, mask='*', dtype='dataset')
        data0 = data[0].T
        data1 = data[1].T

        cpp_data = np.loadtxt('./data/Pr.tz2.dat').transpose()

        aa_eq(data0[0], cpp_data[0]) # distance
        aa_eq(data1[0], cpp_data[0]) # distance
        aa_eq(data0[1], cpp_data[1], decimal=2) # Pr
        aa_eq(data1[1], cpp_data[2], decimal=2) # std

if __name__ == "__main__":
    unittest.main()
