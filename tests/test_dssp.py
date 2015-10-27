#!/usr/bin/env python
import unittest
import numpy as np
import pytraj as pt
from pytraj.testing import aa_eq


class TestDSSP(unittest.TestCase):
    def setUp(self):
        self.traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

    def test_vs_cpptraj(self):
        data = pt.dssp(self.traj, "*", dtype='cpptraj_dataset')
        data_int = np.array([d0.values for d0 in data
                             if d0.dtype == 'integer'],
                            dtype='i4')
        # load cpptraj output
        cpp_data = np.loadtxt("./data/dssp.Tc5b.dat",
                skiprows=1)[:, 1:].T
        aa_eq(data_int.flatten(), cpp_data.flatten())

    def test_frame_indices(self):
        from numpy.testing import assert_equal
        s_0 = pt.dssp(self.traj)[1]
        s_1 = pt.dssp(self.traj, frame_indices=[0, 2, 5])[1]
        assert_equal(s_0[[0, 2, 5]], s_1)

    def test_simplified_codes(self):
        traj = pt.fetch_pdb('1l2y')
        data_full = pt.dssp(traj)[1]
        data_sim = pt.dssp(traj, simplified=True)[1]
        expected_1st = ['C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'C', 'C', 'H', 'H', 'H',
                'H', 'C', 'C', 'C', 'C', 'C', 'C']
        assert expected_1st == data_sim[0].tolist(), 'test_simplified_codes: must equal'

    def test_dssp_all_atoms(self):
        # not assert yet
        traj = pt.fetch_pdb('1l2y')
        data = pt.dssp_all_atoms(traj)

if __name__ == "__main__":
    unittest.main()
