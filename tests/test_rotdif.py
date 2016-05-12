#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.testing import cpptraj_test_dir

cm = """
parm {cpptraj_test_dir}/tz2.parm7
reference {test_rotdif_dir}/avgstruct.pdb [tz2avg] 
trajin {cpptraj_test_dir}/tz2.nc
rms R0 ref [tz2avg] @CA,C,N,O savematrices 
rotdif rmatrix R0[RM] rseed 1 nvecs 10 dt 0.002 tf 0.190 ncorr 101 itmax 500 tol 0.000001 d0 0.03 order 2
""".format(cpptraj_test_dir=cpptraj_test_dir,
           test_rotdif_dir=cpptraj_test_dir + '/Test_Rotdif/')

lines = [line for line in cm.split('\n') if line]

parm = lines[0].split()[-1]
reference = lines[1].split()[1]
trajin = lines[2].split()[-1]

short_cm = ' '.join(lines[-1].split()[1:])
print(short_cm)

class TestRotdif(unittest.TestCase):

    @unittest.skip("cpptraj rotdif does not dump data to DatasetList yet")
    def test_rotdif(self):
        traj = pt.load(trajin, parm)
        ref = pt.load(reference, parm) 

        mat = pt.rotation_matrix(traj, ref=ref, mask='@CA,C,N,O')

        data = pt.all_actions._rotdif(mat, short_cm)

        state = pt.load_cpptraj_state(cm)
        state.run()

if __name__ == "__main__":
    unittest.main()
