#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj import matrix
'''figure out why sign of some eigenvectors are different
'''


class TestDiagMatrix(unittest.TestCase):

    def test_diagmatrix(self):
        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")

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
        indices = np.triu_indices(mat.n_cols)

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
        self.assertRaises(ValueError, lambda: pt.matrix.diagonalize(mat, 3, dtype='ndarray'))

        # plot
        def plot_():
            try:
                from pytraj.plot import plot_matrix
                from matplotlib import pyplot as plt

                fig, ax, bar = plot_matrix(np.abs(data.eigenvectors) - np.abs(
                    cpp_evecs))
                fig.colorbar(bar)
                plt.title('data vs cpp_evecs')

                fig, ax, bar = plot_matrix(np.abs(np_vecs) - np.abs(cpp_evecs))
                fig.colorbar(bar)
                plt.title('np vs cpp_evecs')

                fig, ax, bar = plot_matrix(np.abs(np_vecs) - np.abs(
                    data.eigenvectors))
                fig.colorbar(bar)
                plt.title('np vs data')
                pt.show()
            except ImportError:
                pass


if __name__ == "__main__":
    unittest.main()
