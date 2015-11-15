from __future__ import print_function
import unittest
from collections import OrderedDict
import numpy as np
import pytraj as pt
from pytraj.utils import eq, aa_eq
import pytraj.common_actions as pyca
from pytraj.tools import flatten
from pytraj import matrix
from pytraj.compat import set
from pytraj.parallel import _load_batch_pmap
from pytraj import cpptraj_commands


class TestNormalPmap(unittest.TestCase):
    def setUp(self):
        self.traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

    def test_raise(self):
        # if func is not support pmap
        self.assertRaises(ValueError, lambda: pt.pmap(pt.bfactors, self.traj))

        # run time: openmp vs pmap
        if 'OPENMP' in pt.compiled_info():
            self.assertRaises(RuntimeError,
                              lambda: pt.pmap(pt.watershell, self.traj))

        # if traj is not TrajectoryIterator
        self.assertRaises(ValueError, lambda: pt.pmap(pt.radgyr, self.traj[:]))

        # raise if a given method does not support pmap
        def need_to_raise(traj=self.traj):
            pt.pmap(2, pt.bfactors, traj)

        self.assertRaises(ValueError, lambda: need_to_raise())

        # raise if a traj is not TrajectoryIterator
        def need_to_raise_2(traj=self.traj):
            pt.pmap(pt.bfactors, traj[:], n_cores=2)

        # raise if turn off pmap by setting _is_parallelizable to False
        pt.radgyr._is_parallelizable = False
        self.assertRaises(ValueError, lambda: pt.pmap(pt.radgyr, self.traj))
        pt.radgyr._is_parallelizable = True

    def test_general(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        # with mask
        saved_data = pt.radgyr(traj, '@CA')
        data = pt.pmap(pt.radgyr, traj, '@CA')
        data = pt.tools.dict_to_ndarray(data)
        aa_eq(saved_data, data)

        # with a series of functions
        func_list = [pt.radgyr, pt.molsurf, pt.rmsd]
        ref = traj[-3]

        for n_cores in [2, 3]:
            for func in func_list:
                if func in [pt.rmsd, ]:
                    pout = pt.tools.dict_to_ndarray(pt.pmap(func=func,
                                                            traj=traj,
                                                            ref=ref,
                                                            n_cores=n_cores))
                    serial_out = flatten(func(traj, ref=ref))
                else:
                    pout = pt.tools.dict_to_ndarray(pt.pmap(n_cores=n_cores,
                                                            func=func,
                                                            traj=traj))
                    serial_out = flatten(func(traj))
                aa_eq(pout[0], serial_out)

    def test_different_references(self):
        traj = self.traj
        func = pt.rmsd
        for i in range(0, 8, 2):
            ref = self.traj[i]
            for n_cores in [2, 3, ]:
                pout = pt.tools.dict_to_ndarray(pt.pmap(n_cores=n_cores,
                                                        func=func,
                                                        traj=traj,
                                                        ref=ref))
                serial_out = flatten(func(traj, ref=ref))
                aa_eq(pout[0], serial_out)

    def test_iter_options(self):
        traj = pt.iterload("data/tz2.ortho.nc", "data/tz2.ortho.parm7")
        t0 = traj[:].autoimage().rmsfit(ref=0)
        saved_avg = pt.mean_structure(t0)
        saved_radgyr = pt.radgyr(traj, '@CA')

        # perform autoimage, then rms fit to 1st frame, then compute mean structure
        iter_options = {'autoimage': True, 'rmsfit': 0}
        for n_cores in [2, 3]:
            avg = pt.pmap(pt.mean_structure,
                          traj,
                          iter_options=iter_options,
                          n_cores=n_cores)
            aa_eq(saved_avg.xyz, avg.xyz)
            radgyr_ = pt.tools.dict_to_ndarray(pt.pmap(pt.radgyr,
                                                       traj,
                                                       iter_options={'mask':
                                                                     '@CA'}))
            aa_eq(radgyr_[0], saved_radgyr)


class TestParallelMapForMatrix(unittest.TestCase):
    def test_matrix_module(self):
        traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")

        # not support [covar, distcovar, mwcovar]
        for n_cores in [2, 3]:
            for func in [matrix.dist, matrix.idea]:
                x = pt.pmap(func, traj, '@CA', n_cores=n_cores)
                aa_eq(x, func(traj, '@CA'))

    def test_ired_vector_and_matrix_pmap(self):
        traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")
        h = traj.top.select('@H')
        n = h - 1
        nh = list(zip(n, h))

        exptected_vecs, exptected_mat = pt.ired_vector_and_matrix(traj, nh)
        for n_cores in [2, 3]:
            vecs, mat = pt.pmap(pt.ired_vector_and_matrix,
                                traj,
                                nh,
                                n_cores=n_cores)
            aa_eq(exptected_vecs, vecs, decimal=7)
            aa_eq(exptected_mat, mat, decimal=7)

    def test_rotation_matrix_in_rmsd_calculation(self):
        traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")
        saved_mat = pt.rotation_matrix(traj, ref=traj[3], mask='@CA')
        saved_rmsd = pt.rmsd(traj, ref=traj[3], mask='@CA')

        for n_cores in [2, 3]:
            out = pt.pmap(pt.rotation_matrix, traj, ref=traj[3], mask='@CA')
            out_with_rmsd = pt.pmap(pt.rotation_matrix,
                                    traj,
                                    ref=traj[3],
                                    mask='@CA',
                                    with_rmsd=True)
            mat = out[list(out.keys())[0]]
            mat2, rmsd_ = out_with_rmsd[list(out_with_rmsd.keys())[0]]
            aa_eq(saved_mat, mat)
            aa_eq(saved_mat, mat2)
            aa_eq(saved_rmsd, rmsd_)


class TestCpptrajCommandStyle(unittest.TestCase):
    def test_cpptraj_command_style(self):
        traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")

        angle_ = pt.angle(traj, ':3 :4 :5')
        distance_ = pt.distance(traj, '@10 @20')

        data = pt.pmap(['angle :3 :4 :5', 'distance @10 @20'], traj, n_cores=2)
        assert isinstance(data, OrderedDict), 'must be OrderDict'
        arr = pt.tools.dict_to_ndarray(data)
        aa_eq(angle_, arr[0])
        aa_eq(distance_, arr[1])

        # as whole text, case 1
        data = pt.pmap('''angle :3 :4 :5
        distance @10 @20''',
                       traj,
                       n_cores=2)
        assert isinstance(data, OrderedDict), 'must be OrderDict'
        arr = pt.tools.dict_to_ndarray(data)
        aa_eq(angle_, arr[0])
        aa_eq(distance_, arr[1])

    def test_reference(self):
        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")

        for n_cores in [2, 3]:
            # use 4-th Frame for reference
            data = pt.pmap(['rms @CA refindex 0'],
                           traj,
                           ref=traj[3],
                           n_cores=n_cores)
            arr = pt.tools.dict_to_ndarray(data)[0]
            aa_eq(arr, pt.rmsd(traj, 3, '@CA'))

            # use 4-th and 5-th Frame for reference
            data = pt.pmap(
                ['rms @CA refindex 0', 'rms @CB refindex 1'],
                traj,
                ref=[traj[3], traj[4]],
                n_cores=n_cores)
            arr = pt.tools.dict_to_ndarray(data)[0]
            aa_eq(arr, pt.rmsd(traj, 3, '@CA'))

            arr1 = pt.tools.dict_to_ndarray(data)[1]
            aa_eq(arr1, pt.rmsd(traj, 4, '@CB'))

            # no ref
            data = pt.pmap(['radgyr', ], traj, n_cores=n_cores)
            arr = pt.tools.dict_to_ndarray(data)[0]
            aa_eq(arr, pt.radgyr(traj))


class TestParallelMapForAverageStructure(unittest.TestCase):
    def test_pmap_average_structure(self):
        traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")
        saved_frame = pt.mean_structure(traj, '@CA')
        saved_xyz = saved_frame.xyz

        for n_cores in [2, 3, 4]:
            frame = pt.pmap(pt.mean_structure, traj, '@CA', n_cores=n_cores)
            aa_eq(frame.xyz, saved_xyz)


class TestLoadBathPmap(unittest.TestCase):
    def test_load_batch(self):
        '''just test ValueError
        '''
        self.assertRaises(
            ValueError,
            lambda: _load_batch_pmap(n_cores=4, lines=['autoimage'], traj=None, dtype='dict', root=0, mode='xyz', ref=None))


class TestFrameIndices(unittest.TestCase):
    def test_frame_indices(self):
        traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")

        # frame_indices could be a list, range
        frame_indices_list = [[0, 8, 9, 3, 2, 5], range(6)]

        for frame_indices in frame_indices_list:
            for n_cores in [2, 3]:
                serial_out = pt.radgyr(traj,
                                       '@CA',
                                       frame_indices=frame_indices)
                parallel_out = pt.pmap(pt.radgyr,
                                       traj,
                                       '@CA',
                                       frame_indices=frame_indices)
                parallel_out_cpptraj_style = pt.pmap(
                    ['radgyr @CA nomax'],
                    traj,
                    frame_indices=frame_indices)
                aa_eq(serial_out, pt.tools.dict_to_ndarray(parallel_out))
                aa_eq(serial_out,
                      pt.tools.dict_to_ndarray(parallel_out_cpptraj_style))


class TestCheckValidCommand(unittest.TestCase):
    def test_check_valid_command(self):
        traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")

        # must provide refindex
        self.assertRaises(ValueError, lambda: pt.pmap(['rms'], traj, n_cores=2))
        self.assertRaises(ValueError, lambda: pt.pmap(('rms',), traj, n_cores=2))
        self.assertRaises(ValueError, lambda: pt.pmap('rms', traj, n_cores=2))
        # does not support matrix
        self.assertRaises(ValueError, lambda: pt.pmap(['matrix'], traj, n_cores=2))

        # do not accept any cpptraj analysis command
        for word in cpptraj_commands.analysis_commands:
            self.assertRaises(ValueError, lambda: pt.pmap(word, traj, n_cores=2))


if __name__ == "__main__":
    unittest.main()
