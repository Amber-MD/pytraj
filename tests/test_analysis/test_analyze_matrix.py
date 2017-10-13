#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import unittest
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq
from pytraj import matrix
'''figure out why sign of some eigenvectors are different
'''


class TestDiagMatrix(unittest.TestCase):
    def test_diagmatrix(self):
        traj = pt.iterload(fn('tz2.nc'), fn('tz2.parm7'))

        state = pt.load_batch(traj, '''
        #matrix covar @CA name mymat
        matrix dist @CA name mymat
        diagmatrix mymat vecs 6''')

        state.run()
        mat = state.data['mymat']
        cpp_evals, cpp_evecs = state.data[-1].eigenvalues, state.data[
            -1].eigenvectors

        # test triu_indices
        mat2 = mat.__class__()
        np.triu_indices(mat.n_cols)

        mat2._set_data_half_matrix(mat._to_cpptraj_sparse_matrix(), mat.size,
                                   mat.n_cols)

        # OK
        data = matrix.diagonalize(mat.values, n_vecs=6, dtype='dataset')[-1]
        data2 = matrix.diagonalize(mat, n_vecs=6, dtype='dataset')[0]
        aa_eq(data.eigenvalues, cpp_evals, decimal=10)
        aa_eq(data2.eigenvalues, cpp_evals, decimal=10)

        # try numpy
        np_vals, np_vecs = np.linalg.eigh(mat2.values, UPLO='U')
        np_vals = np_vals[::-1][:6]
        np_vecs = np_vecs[:, ::-1].T[:6]

        # at least the absolute values are equal
        aa_eq(np.abs(data.eigenvectors), np.abs(cpp_evecs))
        aa_eq(np.abs(np_vecs), np.abs(cpp_evecs))

        # test raise if not having supported dtype
        self.assertRaises(
            ValueError, lambda: pt.matrix.diagonalize(mat, 3, dtype='ndarray'))

    def test_diagmatrix_mwcovar(self):
        traj = pt.iterload(fn('tz2.nc'), fn('tz2.parm7'))

        state = pt.load_batch(traj, '''
        matrix mwcovar @CA name mymat
        diagmatrix mymat vecs 6 name mydiag''')
        state.run()

        mat = pt.matrix.mwcovar(traj, '@CA')
        ca_indices = traj.top.select('@CA')
        eigenvectors, eigenvalues = pt.matrix.diagonalize(
            mat,
            n_vecs=6,
            scalar_type='mwcovar',
            mass=traj.top.mass[ca_indices])
        aa_eq(np.abs(state.data['mydiag'].eigenvalues), np.abs(eigenvalues))
        aa_eq(np.abs(state.data['mydiag'].eigenvectors), np.abs(eigenvectors))


if __name__ == "__main__":
    unittest.main()
