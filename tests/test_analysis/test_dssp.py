#!/usr/bin/env python
import unittest
import numpy as np
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq

from utils import fn


def print_name(func):
    print(func.__name__)
    return func


class TestDSSP(unittest.TestCase):
    @print_name
    def test_vs_cpptraj(self):
        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))
        data = pt.dssp(traj, "*", dtype='cpptraj_dataset')
        data_int = np.array(
            [d0.values for d0 in data if d0.dtype == 'integer'], dtype='i4')
        # load cpptraj output
        cpp_data = np.loadtxt(fn("dssp.Tc5b.dat"), skiprows=1)[:, 1:].T
        aa_eq(data_int.flatten(), cpp_data.flatten())

    @print_name
    def test_frame_indices(self):
        from numpy.testing import assert_equal
        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))
        s_0 = pt.dssp(traj)[1]
        s_1 = pt.dssp(traj, frame_indices=[0, 2, 5])[1]
        assert_equal(s_0[[0, 2, 5]], s_1)

    @print_name
    def test_simplified_codes(self):
        traj = pt.load(fn("1L2Y.pdb"))
        pt.dssp(traj)[1]
        data_sim = pt.dssp(traj, simplified=True)[1]
        expected_1st = [
            'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'C', 'C', 'H', 'H', 'H',
            'H', 'C', 'C', 'C', 'C', 'C', 'C'
        ]
        assert expected_1st == data_sim[
            0].tolist(), 'test_simplified_codes: must equal'

    @print_name
    def test_dssp_allatoms(self):
        traj = pt.load(fn("1L2Y.pdb"))
        allatoms_dssp = pt.dssp_allatoms(traj).T
        allresidues_dssp = pt.dssp_allresidues(traj).T
        for idx, res in enumerate(traj.top.residues):
            assert np.all(
                np.array(allatoms_dssp[res.first_atom_index:
                                       res.last_atom_index]) == (
                                           np.array(allresidues_dssp[idx])))


if __name__ == "__main__":
    unittest.main()
