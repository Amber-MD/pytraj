#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq

command = '''
# Step one. Generate average structure.
# RMS-Fit to first frame to remove global translation/rotation.
parm tz2.parm7
trajin tz2.nc
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


class TestProject(unittest.TestCase):

    def test_projection(self):
        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")


if __name__ == "__main__":
    unittest.main()
