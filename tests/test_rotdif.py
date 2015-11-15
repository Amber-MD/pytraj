#!/usr/bin/env python
from __future__ import print_function
import os
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.testing import cpptraj_test_dir


class TestRotdif(unittest.TestCase):
    # TODO: make sure cpptraj saving results to Dataset
    # not assertion yet

    @unittest.skip('skip rotdif, too verbose')
    def test_rotdif(self):
        '''test_rotdif
        '''
        avg_fn = os.path.join(cpptraj_test_dir, 'Test_Rotdif', 'avgstruct.pdb')
        traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")
        avg = pt.iterload(avg_fn, traj.top)
        mat = pt.calc_rotation_matrix(traj, ref=avg, mask='@CA,C,N,O')
        data = pt._rotdif(mat,
                          rseed=1,
                          nvecs=10,
                          dt=0.002,
                          tf=0.19,
                          itmax=500,
                          tol=0.000001,
                          d0=0.03,
                          order=2)

        text = '''
        parm data/tz2.parm7
        reference %s [tz2avg]
        trajin data/tz2.nc
        rms R0 ref [tz2avg] @CA,C,N,O savematrices
        rotdif rmatrix R0[RM] rseed 1 nvecs 10 dt 0.002 tf 0.190        itmax 500 tol 0.000001 d0
        0.03 order 2''' % avg_fn
        # too verbose


if __name__ == "__main__":
    unittest.main()
