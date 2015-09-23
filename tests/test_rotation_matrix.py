#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")

        # cpptraj output
        cm = '''
        average avg.pdb
        reference avg.pdb
        rms R0 reference @CA,C,N,O savematrices
        '''
        s = pt.load_batch(traj, cm)
        s.run()
        saved_mat = s.data[-1].values

        # pytraj output
        avg = pt.mean_structure(traj)
        mat = pt.calc_rotation_matrix(traj, ref=avg, mask='@CA,C,N,O')
        assert mat.shape == (traj.n_frames, 3, 3), 'mat shape'
        aa_eq(mat.flatten(), saved_mat.flatten())


if __name__ == "__main__":
    unittest.main()
