#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestNHOrderParamters(unittest.TestCase):
    def test_nh_paramters(self):
        parmfile =  '../cpptraj/test/Test_IRED/1IEE_A_prot.prmtop'
        trajfile = '../cpptraj/test/Test_IRED/1IEE_A_test.mdcrd'
        traj = pt.iterload(trajfile, parmfile)

        h_indices = pt.select_atoms(traj.top, '@H')
        n_indices = h_indices - 1
        nh_indices = list(zip(n_indices, h_indices))

        orders = pt.NH_order_paramters(traj, nh_indices)
        saved_S2 = np.loadtxt('../cpptraj/test/Test_IRED/orderparam.save').T[-1]

        aa_eq(orders, saved_S2)


if __name__ == "__main__":
    unittest.main()
