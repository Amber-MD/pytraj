#!/usr/bin/env python
import unittest
import numpy as np
import pytraj as pt
from pytraj.testing import aa_eq


class TestDSSP(unittest.TestCase):
    def test_vs_cpptraj(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        data = pt.dssp(traj, "*", dtype='cpptraj_dataset')
        data_int = np.array([d0.values for d0 in data
                             if d0.dtype == 'integer'],
                            dtype='i4')
        # load cpptraj output
        cpp_data = np.loadtxt("./data/dssp.Tc5b.dat",
                              skiprows=1).transpose()[1:]
        aa_eq(data_int.flatten(), cpp_data.flatten())

    def test_frame_indices(self):
        from numpy.testing import assert_equal
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        s_0 = pt.dssp(traj)[1]
        s_1 = pt.dssp(traj, frame_indices=[0, 2, 5])[1]
        assert_equal(s_0[:, [0, 2, 5]], s_1)


if __name__ == "__main__":
    unittest.main()
