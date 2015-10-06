#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj import matrix


class TestDiagMatrix(unittest.TestCase):
    def test_diagmatrix(self):
        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")
        state = pt.load_batch(traj, '''
        matrix covar name mat
        diagmatrix mat vecs 6''')
        state.run()
        mat = state.data[1]
        evals, evecs = state.data[-1].eigenvalues, state.data[-1].eigenvectors
        #print(mat, mat.kind)
        #print(evecs, evecs)

        data = matrix.diagonalize(mat, nvecs=6)
        #print(data)


if __name__ == "__main__":
    unittest.main()
