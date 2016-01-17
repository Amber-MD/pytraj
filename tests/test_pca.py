#!/usr/bin/env python

from __future__ import print_function
import unittest
import numpy as np
import pytraj as pt
from pytraj.utils import eq, aa_eq

command = '''
# Step one. Generate average structure.
# RMS-Fit to first frame to remove global translation/rotation.
parm data/tz2.parm7
trajin data/tz2.nc
rms first !@H=
average crdset AVG
run
# Step two. RMS-Fit to average structure. Calculate covariance matrix.
# Save the fit coordinates.
rms ref AVG !@H=
matrix covar name MyMatrix !@H=
createcrd CRD1
run
# Step three. Diagonalize matrix.
runanalysis diagmatrix MyMatrix vecs 2 name MyEvecs
# Step four. Project saved fit coordinates along eigenvectors 1 and 2
crdaction CRD1 projection evecs MyEvecs !@H= out project.dat beg 1 end 2
'''

command_ref_provided = '''
parm data/tz2.parm7
reference data/tz2.rst7
trajin data/tz2.nc
rms reference !@H=
matrix covar name MyMatrix !@H=
createcrd CRD1
run
# Step three. Diagonalize matrix.
runanalysis diagmatrix MyMatrix vecs 2 name MyEvecs
# Step four. Project saved fit coordinates along eigenvectors 1 and 2
crdaction CRD1 projection evecs MyEvecs !@H= out project.dat beg 1 end 2
'''


class TestPCA(unittest.TestCase):

    def test_pca_noref(self):
        '''no reference'''
        traj = pt.load("data/tz2.nc", "data/tz2.parm7")

        # no reference
        state = pt.load_cpptraj_state(command)
        state.run()

        mask = '!@H='

        data = pt.pca(traj, mask, n_vecs=2)
        cpp_data = state.data[-2:].values
        # use absolute values
        aa_eq(np.abs(data[0]), np.abs(cpp_data), decimal=3)

    def test_pca_with_ref(self):
        '''has reference'''
        traj = pt.load("data/tz2.nc", "data/tz2.parm7")
        ref = pt.load('data/tz2.rst7', traj.top)

        state = pt.load_cpptraj_state(command_ref_provided)
        state.run()

        mask = '!@H='

        data = pt.pca(traj, mask, n_vecs=2, ref=ref)
        cpp_data = state.data[-2:].values
        # use absolute values
        aa_eq(np.abs(data[0]), np.abs(cpp_data), decimal=3)

    def test_pca_raise(self):
        traj = pt.iterload('data/tz2.nc', 'data/tz2.parm7')
        self.assertRaises(ValueError, lambda: pt.pca(traj, n_vecs=2, mask='@CA'))


if __name__ == "__main__":
    unittest.main()
