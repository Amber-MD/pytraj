#!/usr/bin/env python

from __future__ import print_function
import unittest
import numpy as np
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestPCA(unittest.TestCase):

    def test_pca_noref(self):
        '''no reference
        
        pytraj: pt.pca(traj, mask, n_vecs=2) 
        '''

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

        traj = pt.load("data/tz2.nc", "data/tz2.parm7")

        # no reference
        state = pt.load_cpptraj_state(command)
        state.run()

        mask = '!@H='

        data = pt.pca(traj, mask, n_vecs=2)
        cpp_data = state.data[-2:].values
        # use absolute values
        aa_eq(np.abs(data[0]), np.abs(cpp_data), decimal=3)

    def test_pca_noref_nofit(self):
        '''no reference and do not do fitting
        
        from drroe: " Also, not fitting at all should be considered a legitimate option - 
        you may want to include global rotational and translational motion in your eigenvectors."

        pytraj: pt.pca(traj, mask, n_vecs=2, fit=False)
        '''

        command = '''
        parm data/tz2.parm7
        trajin data/tz2.nc

        matrix covar name MyMatrix !@H=
        createcrd CRD1
        run

        # Step three. Diagonalize matrix.
        runanalysis diagmatrix MyMatrix vecs 2 name MyEvecs

        # Step four. Project saved fit coordinates along eigenvectors 1 and 2
        crdaction CRD1 projection evecs MyEvecs !@H= out project.dat beg 1 end 2
        '''

        traj = pt.load("data/tz2.nc", "data/tz2.parm7")

        # no reference
        state = pt.load_cpptraj_state(command)
        state.run()

        mask = '!@H='

        data = pt.pca(traj, mask, n_vecs=2, fit=False)
        data_ref = pt.pca(traj, mask, n_vecs=2, fit=False, ref=3, ref_mask='@CA')
        cpp_data = state.data[-2:].values
        # use absolute values
        aa_eq(np.abs(data[0]), np.abs(cpp_data), decimal=3)
        # if fit=True, ref will be ignored
        aa_eq(np.abs(data_ref[0]), np.abs(cpp_data), decimal=3)

    def test_pca_with_ref(self):
        '''has reference

        from drroe: "If the user provides their own reference structure, do not create an average structure"

        pytraj: pt.pca(traj, mask, n_vecs=2, ref=ref)
        '''

        command_ref_provided = '''
        parm data/tz2.parm7
        trajin data/tz2.nc
        reference data/tz2.rst7

        rms reference !@H=

        matrix covar name MyMatrix !@H=
        createcrd CRD1
        run

        # Step three. Diagonalize matrix.
        runanalysis diagmatrix MyMatrix vecs 2 name MyEvecs

        # Step four. Project saved fit coordinates along eigenvectors 1 and 2
        crdaction CRD1 projection evecs MyEvecs !@H= out project.dat beg 1 end 2
        '''

        traj = pt.load("data/tz2.nc", "data/tz2.parm7")
        ref = pt.load('data/tz2.rst7', traj.top)

        state = pt.load_cpptraj_state(command_ref_provided)
        state.run()

        mask = '!@H='

        data = pt.pca(traj, mask, n_vecs=2, ref=ref)
        cpp_data = state.data[-2:].values

        # use absolute values
        aa_eq(np.abs(data[0]), np.abs(cpp_data), decimal=3)

    def test_pca_with_ref_with_different_mask_from_matrix(self):
        '''has reference. Use !@H= for ref_mask and use * for covariance matrix  and projection

        from drroe: "You should be able to supply separate masks for fitting and creating the covariance matrix
        It is common enough for example to only perform rms-fitting on heavy atoms while still wanting all atoms in eigenvectors."

        pytraj: pt.pca(traj, mask=mask_matrix, n_vecs=2, ref=ref, ref_mask=mask_ref)
        '''

        command_ref_provided = '''
        parm data/tz2.parm7
        trajin data/tz2.nc

        reference data/tz2.rst7

        # only perform fitting on heavy atoms
        rms reference !@H=

        # all atoms
        matrix covar name MyMatrix *
        createcrd CRD1
        run

        # Step three. Diagonalize matrix.
        runanalysis diagmatrix MyMatrix vecs 2 name MyEvecs

        # Step four. Project saved fit coordinates along eigenvectors 1 and 2
        # all atoms
        crdaction CRD1 projection evecs MyEvecs * out project.dat beg 1 end 2
        '''

        traj = pt.load("data/tz2.nc", "data/tz2.parm7")
        ref = pt.load('data/tz2.rst7', traj.top)

        state = pt.load_cpptraj_state(command_ref_provided)
        state.run()

        mask_ref = '!@H='
        mask_matrix = '*'

        data = pt.pca(traj, mask=mask_matrix, n_vecs=2, ref=ref, ref_mask=mask_ref)
        cpp_data = state.data[-2:].values
        # use absolute values
        aa_eq(np.abs(data[0]), np.abs(cpp_data), decimal=3)

    def test_pca_raise(self):
        traj = pt.iterload('data/tz2.nc', 'data/tz2.parm7')
        self.assertRaises(ValueError, lambda: pt.pca(traj, n_vecs=2, mask='@CA'))


if __name__ == "__main__":
    unittest.main()
