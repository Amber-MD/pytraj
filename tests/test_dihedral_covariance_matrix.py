#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestDihedralCovarianceMatrix(unittest.TestCase):
    # TODO: add assertion

    def test_dihcovar(self):
        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")
        txt = '''
        parm data/tz2.parm7
        trajin data/tz2.nc

        # Generation of phi/psi dihedral data
        multidihedral BB phi psi resrange 2

        run
        # Calculate dihedral covariance matrix and obtain eigenvectors
        matrix dihcovar dihedrals BB[*] name DIH

        diagmatrix DIH vecs 4 name DIHMODES
        run
        # Project along eigenvectors
        projection evecs DIHMODES beg 1 end 4 dihedrals BB[*]
        '''

        state = pt.load_cpptraj_state(txt)
        state.run()

        modes = state.data['DIHMODES']
        # need to provide all dihedrals too?


if __name__ == "__main__":
    unittest.main()
