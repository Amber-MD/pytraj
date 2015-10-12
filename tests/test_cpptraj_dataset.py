#!/usr/bin/env python

from __future__ import print_function
import unittest
import numpy as np
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.datasets import cpp_datasets
from pytraj.datasets import CpptrajDatasetList

# used for loading cpptraj state
txt = '''
# Step one. Generate average structure.
# RMS-Fit to first frame to remove global translation/rotation.
parm data/tz2.parm7
trajin data/tz2.nc
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

#@unittest.skipIf('DNO_MATHLIB' in pt.compiled_info(), 'there is no LAPACK')
@unittest.skip('just skip')
class TestCpptrajDatasetWithMathLib(unittest.TestCase):
    def setUp(self):
        self.state = pt.datafiles.load_cpptraj_state(txt)
        self.state.run()
        self.traj = pt.iterload('data/tz2.nc', 'data/tz2.parm7')

    def test_call_values(self):
        for d in self.state.data:
            d.values

    def test_dataset_coords_ref(self):
        traj = pt.iterload('data/tz2.nc', 'data/tz2.parm7')
        avg_frame = pt.mean_structure(traj(rmsfit=(0, '!@H=')))
        state = self.state

        # need to loop several times to make sure this does not fail
        # due to memory free
        for _ in range(20):
            cpp_ref = state.data['AVG'].get_frame()
            aa_eq(avg_frame.xyz, cpp_ref.xyz)

    def test_DatasetModes(self):
        state = self.state
        modes = state.data['MyEvecs']
        mat = state.data['MyMatrix'].values
        np_eg = np.linalg.eigh(mat)

        # make sure eigenvalues from cpptraj are the same as ones in numpy
        # we calculated only 2 eigenvalues
        aa_eq(sorted(modes.eigenvalues), np_eg[0][-2:])
        aa_eq(modes.eigenvectors[0], np_eg[1][:, -1])
        aa_eq(modes.eigenvectors[1], np_eg[1][:, -2])

class TestCpptrajDatasetWithoutMathLib(unittest.TestCase):
    def setUp(self):
        self.traj = pt.iterload('data/tz2.nc', 'data/tz2.parm7')

    def test_DatasetMatrix3x3(self):
        # test _append_from_array
        mat0 = pt.calc_rotation_matrix(self.traj, ref=0)

        shape2d = (mat0.shape[0], mat0.shape[1] * mat0.shape[2])
        dmat3x3 = cpp_datasets.DatasetMatrix3x3()
        dmat3x3._append_from_array(mat0.reshape(shape2d))
        aa_eq(mat0, dmat3x3.values)

    def test_DatasetMatrixDouble(self):
        from pytraj import matrix

        for op in [matrix.dist, matrix.covar]:
            orig_mat = op(self.traj, '!@H=', dtype='cpptraj_dataset')[0]
            shape = orig_mat.values.shape
            cpp_mat = orig_mat._to_cpptraj_sparse_matrix()
            assert orig_mat.kind == 'half', 'must be half matrix'

            new_mat = orig_mat.__class__()
            new_mat._set_data_half_matrix(cpp_mat, orig_mat.size, shape[0])
            assert new_mat.kind == 'half', 'new_mat must be half matrix'
            aa_eq(orig_mat.values, new_mat.values)

    def test_add_new_for_CpptrajDatasetList(self):
        #"MATRIX_DOUBLE" : MATRIX_DBL, 
        #"MATRIX_FLOAT" : MATRIX_FLT,
        #"VECTOR" : VECTOR,
        #"MODES" : MODES,
        #"GRID_FLOAT" : GRID_FLT,
        #"XYMESH" : XYMESH, 
        #"TRAJ" : TRAJ, 
        dslist = CpptrajDatasetList()
        
        # integer
        dslist.add_new('integer', name='my_int')
        dslist[-1].data = [2, 3]
        aa_eq(dslist[-1].values, [2, 3])

        # double
        dslist.add_new('double', name='my_double')
        dslist[-1].data = [2, 3]
        aa_eq(dslist[-1].values, [2, 3])

        # double
        dslist.add_new('float', name='my_float')
        dslist[-1].data = [2, 3]
        aa_eq(dslist[-1].values, [2, 3])

        # string
        dslist.add_new('string', name='my_string')
        dslist[-1].data = ['H', 'T']
        assert dslist[-1].values.tolist() == ['H', 'T'], 'string must be equal'

        # reference
        dslist.add_new('reference', name='my_reference')
        dslist[-1].data = self.traj[-2]
        aa_eq(dslist[-1].xyz, self.traj[-2].xyz)

        # matrix3x3
        dslist.add_new('matrix3x3', name='my_mat3x3')
        mat = pt.calc_rotation_matrix(self.traj, ref=0, mask='@CA')
        # there is no assignment. Need to update by another method
        dslist[-1]._append_from_array(mat)
        aa_eq(dslist[-1].values, mat)

        # TRAJ
        dslist.add_new('traj', name='my_traj')
        dslist[-1].top = self.traj.top
        dslist[-1].load(self.traj.filename)
        traj_new = dslist[-1]
        traj_new.top._own_memory = False
        print(traj_new._own_memory, traj_new.top._own_memory)
        print(traj_new[::2])
        #for _ in traj_new: pass
        #aa_eq(traj_new.xyz, self.traj.xyz)


if __name__ == "__main__":
    unittest.main()
