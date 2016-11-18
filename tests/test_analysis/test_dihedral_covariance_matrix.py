#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from utils import fn, tz2_trajin, tz2_top


class TestDihedralCovarianceMatrix(unittest.TestCase):
    # TODO: add assertion

    def test_dihcovar(self):
        pt.iterload(fn('tz2.nc'), fn('tz2.parm7'))
        txt = '''
        parm {}
        trajin {}

        # Generation of phi/psi dihedral data
        multidihedral BB phi psi resrange 2

        run
        # Calculate dihedral covariance matrix and obtain eigenvectors
        matrix dihcovar dihedrals BB[*] name DIH

        diagmatrix DIH vecs 4 name DIHMODES
        run
        # Project along eigenvectors
        projection evecs DIHMODES beg 1 end 4 dihedrals BB[*]
        '''.format(tz2_top, tz2_trajin)

        state = pt.load_cpptraj_state(txt)
        state.run()

        state.data['DIHMODES']
        # need to provide all dihedrals too?


if __name__ == "__main__":
    unittest.main()
