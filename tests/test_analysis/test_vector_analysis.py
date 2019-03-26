#!/usr/bin/env python

from __future__ import print_function
import unittest
import os
import numpy as np
import pytraj as pt
from utils import fn
from pytraj import vector as va
from pytraj.testing import aa_eq
from pytraj.testing import cpptraj_test_dir
from pytraj.datasets.c_datasetlist import DatasetList
from pytraj import ActionList
from pytraj.analysis.c_action import c_action as CA


class TestVectorAnalysisModule(unittest.TestCase):
    def test_vector_principal(self):
        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))
        pdata = pt.vector.principal(traj)
        assert pdata.shape == (10, 3)
        np.testing.assert_almost_equal(
            pt.vector.principal(traj),
            [[-0.16903813783287722, 0.36316033000258224, 0.9162645265808389],
             [-0.3199676454090604, 0.8725925204063159, -0.36905690512756323], [
                 0.9845639960055733, -0.1655283170881773, -0.056869271241093346
             ], [0.806423302456264, 0.35877569174107243, -0.47006537872007553],
             [0.06177153699316464, 0.9689358562119692, 0.23947355545922028],
             [-0.1977750666117889, 0.10362572565059913, 0.9747547035075487],
             [0.8122317492961045, -0.032170512063985766, 0.5824471165598908], [
                 0.30670002856835077, 0.22464432292840247, -0.9249162235857972
             ], [0.7710019864344335, -0.04553437498669517, 0.6352027688925233],
             [-0.17181943718191456, 0.7403546691757756, 0.6498869477348488]])

    def test_vector_raise(self):
        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))
        arr = np.arange(100).astype('i4').reshape(2, 25, 2)
        self.assertRaises(ValueError, lambda: pt.vector.vector_mask(traj, arr))

    def test_actionlist(self):
        '''test_actionlist
        '''
        dslist = DatasetList()
        actlist = ActionList()
        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))
        mask_list = ['@CB @CA', '@CA @H']

        for mask in mask_list:
            actlist.add(CA.Action_Vector(), mask, traj.top, dslist=dslist)
        actlist.compute(traj)

        dslist4 = va.vector_mask(traj, mask_list)

        dslist3_0 = va.vector_mask(traj, mask_list[0])
        dslist3_1 = va.vector_mask(traj, mask_list[1])

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

        with open(fn('ired.in'), 'r') as fh:
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
        cpp_data = state.data
        cpp_vectors = cpp_data.grep('vector', mode='dtype').values
        cpp_matired = cpp_data.grep('matrix', mode='dtype')['matired']

        # assert between pytraj's data_vec and cpptraj's cpp_vectors
        aa_eq(data_vec, cpp_vectors)

        # from a 2D array of integers
        aa_eq(data_vec_2, cpp_vectors)

        # test ired vector with ired matrix
        # open file

        with open(fn('ired_reduced.in'), 'r') as fh:
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
