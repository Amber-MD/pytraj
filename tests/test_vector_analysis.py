#!/usr/bin/env python

from __future__ import print_function
import unittest
import os
import numpy as np
import pytraj as pt
from pytraj import analdict
from pytraj import vector as va
from pytraj import adict
from pytraj.testing import aa_eq
from pytraj.testing import cpptraj_test_dir
from pytraj.datasets.c_datasetlist import DatasetList
from pytraj import ActionList
from pytraj.c_action import c_action as CA


class TestVectorAnalysisModule(unittest.TestCase):

    def test_vector_raise(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        arr = np.arange(100).astype('i4').reshape(2, 25, 2)
        self.assertRaises(ValueError, lambda: pt.vector.vector_mask(traj, arr))

    def test_actionlist(self):
        '''test_actionlist
        '''
        dslist = DatasetList()
        actlist = ActionList()
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        mask_list = ['@CB @CA', '@CA @H']

        for mask in mask_list:
            actlist.add(CA.Action_Vector(),
                        mask,
                        traj.top,
                        dslist=dslist)
        actlist.compute(traj)

        dslist2 = pt.calc_vector(traj, mask_list)
        dslist4 = va.vector_mask(traj, mask_list)

        dslist3_0 = pt.calc_vector(traj, mask_list[0])
        dslist3_1 = pt.calc_vector(traj, mask_list[1])

        aa_eq(dslist3_0, dslist2[0])
        aa_eq(dslist3_1, dslist2[1])

        aa_eq(dslist3_0, dslist4[0])
        aa_eq(dslist3_1, dslist4[1])

        aa_eq(dslist3_0, dslist[0].values)
        aa_eq(dslist3_1, dslist[1].values)

    def test_ired_vector(self):
        '''test mask as a list of strings or as a 2D array of integers
        '''
        parm_dir = os.path.join(cpptraj_test_dir, 'Test_IRED',
                                '1IEE_A_prot.prmtop')
        trajin_dir = os.path.join(cpptraj_test_dir, 'Test_IRED',
                                  '1IEE_A_test.mdcrd')
        traj = pt.iterload(trajin_dir, parm_dir)

        # get a list of mask from cpptraj input
        maskes = []
        lines = None

        n_indices_cpp = []

        with open('data/ired.in', 'r') as fh:
            lines = fh.readlines()
            for line in lines:
                if 'vector' in line and 'ired' in line:
                    # example: vector v100 @1541 ired @1542
                    sline = line.split()
                    mask = ' '.join((sline[2], sline[4]))
                    n_indices_cpp.append(int(sline[2][1:]) - 1)
                    maskes.append(mask)

        h_indices_cpp = [i + 1 for i in n_indices_cpp]

        # calcuate vector from a list of strings
        data_vec = va.vector_mask(traj, maskes)

        # calcuate vector from a 2d array of integers
        nh_indices = np.array(list(zip(n_indices_cpp, h_indices_cpp)))
        data_vec_2 = va.vector_mask(traj, nh_indices)

        # re-create cpptraj input to run cpptraj
        txt = ''.join(lines)
        # add parm and trajin lines
        txt = 'parm ' + parm_dir + '\n' + \
              'trajin ' + trajin_dir + '\n' + \
              txt

        state = pt.datafiles.load_cpptraj_output(txt, dtype='state')
        state.run()
        cpp_data = state.datasetlist
        cpp_vectors = cpp_data.grep('vector', mode='dtype').values
        cpp_matired = cpp_data.grep('matrix', mode='dtype')['matired']

        # assert between pytraj's data_vec and cpptraj's cpp_vectors
        aa_eq(data_vec, cpp_vectors)

        # from a 2D array of integers
        aa_eq(data_vec_2, cpp_vectors)

        # test ired vector with ired matrix
        # open file

        with open('data/ired_reduced.in', 'r') as fh:
            text = ''.join(fh.readlines())
        state2 = pt.load_batch(traj, text)
        state2.run()

        data = pt.ired_vector_and_matrix(traj, nh_indices, order=2)
        data_vec_3 = data[0]
        assert len(data_vec_3) == 126, 'must have 126 vectors'
        matired = data[1]
        # TODO: know why??
        matired /= matired[0, 0]
        aa_eq(data_vec_3, cpp_vectors)
        assert pt.tools.rmsd(matired.flatten(),
                             cpp_matired.values) < 1E-6, 'matired'


if __name__ == "__main__":
    unittest.main()
