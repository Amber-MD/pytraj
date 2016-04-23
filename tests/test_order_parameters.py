#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import sys
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj import pmap
from pytraj.testing import cpptraj_test_dir


class TestNHOrderParamters(unittest.TestCase):

    @unittest.skipIf(sys.platform != 'linux', 'pmap for linux')
    def test_nh_paramters(self):
        parmfile = cpptraj_test_dir + '/Test_IRED/1IEE_A_prot.prmtop'
        trajfile = cpptraj_test_dir + '/Test_IRED/1IEE_A_test.mdcrd'
        traj = pt.iterload(trajfile, parmfile)

        h_indices = pt.select_atoms('@H', traj.top)
        n_indices = h_indices - 1
        nh_indices = list(zip(n_indices, h_indices))

        # single core
        orders = pt.NH_order_parameters(traj, nh_indices, tcorr=8000.)
        saved_S2 = np.loadtxt(cpptraj_test_dir +
                              '/Test_IRED/orderparam.save').T[-1]

        aa_eq(orders, saved_S2)

        # multiple core
        # default 2
        orders = pt.pmap(pt.NH_order_parameters, traj, nh_indices, tcorr=8000.)
        aa_eq(orders, saved_S2)

        for n_cores in [1, 2, 3, 4, -1]:
            orders = pt.NH_order_parameters(traj,
                                            nh_indices,
                                            tcorr=8000.,
                                            n_cores=n_cores)
            aa_eq(orders, saved_S2)

            orders = pt.pmap(pt.NH_order_parameters,
                             traj,
                             nh_indices,
                             tcorr=8000.,
                             n_cores=n_cores)
            aa_eq(orders, saved_S2)


if __name__ == "__main__":
    unittest.main()
