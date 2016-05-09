#!/usr/bin/env python
import unittest
import numpy as np
import pytraj as pt
from pytraj.testing import aa_eq

try:
    import mdtraj as md
    has_mdtraj = True
except ImportError:
    has_mdtraj = False


class TestDSSP(unittest.TestCase):

    def setUp(self):
        self.traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")

    def test_vs_cpptraj(self):
        data = pt.dssp(self.traj, "*", dtype='cpptraj_dataset')
        data_int = np.array(
            [d0.values for d0 in data if d0.dtype == 'integer'],
            dtype='i4')
        # load cpptraj output
        cpp_data = np.loadtxt("./data/dssp.Tc5b.dat", skiprows=1)[:, 1:].T
        aa_eq(data_int.flatten(), cpp_data.flatten())

    def test_frame_indices(self):
        from numpy.testing import assert_equal
        s_0 = pt.dssp(self.traj)[1]
        s_1 = pt.dssp(self.traj, frame_indices=[0, 2, 5])[1]
        assert_equal(s_0[[0, 2, 5]], s_1)

    def test_simplified_codes(self):
        traj = pt.load("data/1L2Y.pdb")
        data_full = pt.dssp(traj)[1]
        data_sim = pt.dssp(traj, simplified=True)[1]
        expected_1st = ['C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'C', 'C', 'H',
                        'H', 'H', 'H', 'C', 'C', 'C', 'C', 'C', 'C']
        assert expected_1st == data_sim[0].tolist(
        ), 'test_simplified_codes: must equal'

    @unittest.skipIf(not has_mdtraj, 'need mdtraj to assert')
    def test_dssp_allresidues(self):
        from numpy.testing import assert_array_equal

        def update_mdtraj_dssp(mdata):
            for idx, elm in enumerate(mdata):
                if elm == 'NA':
                    mdata[idx] = 'C'
            return mdata

        trajlist = []
        trajlist.append(pt.iterload('data/DPDP.nc', 'data/DPDP.parm7'))
        trajlist.append(pt.iterload('data/tz2.ortho.nc',
                                    'data/tz2.ortho.parm7'))

        trajlist.append(pt.iterload('data/1L2Y.pdb'))

        for traj in trajlist:
            data = pt.dssp_allresidues(traj, simplified=True)[0]

            mtraj = md.load(traj.filename, top=traj.top.filename)
            mdata = md.compute_dssp(mtraj, simplified=True)[0]
            mdata = update_mdtraj_dssp(mdata)
            assert_array_equal(data, mdata)

    def test_dssp_allatoms(self):
        traj = pt.load("data/1L2Y.pdb")
        allatoms_dssp = pt.dssp_allatoms(traj).T
        allresidues_dssp = pt.dssp_allresidues(traj).T
        for idx, res in enumerate(traj.top.residues):
            assert np.all(np.array(allatoms_dssp[res.first_atom_index: res.last_atom_index])
                          == (np.array(allresidues_dssp[idx])))

if __name__ == "__main__":
    unittest.main()
