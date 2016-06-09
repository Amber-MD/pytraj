from __future__ import print_function
import sys
import unittest
from collections import OrderedDict
import numpy as np
import pytraj as pt
from pytraj.utils import eq, aa_eq

from pytraj.tools import flatten
from pytraj import matrix
from pytraj.compat import set
from pytraj.parallel.base import _load_batch_pmap, worker_by_actlist
from pytraj.parallel.multiprocess import worker_byfunc
from pytraj import c_commands


@unittest.skipUnless(sys.platform.startswith('linux'), 'pmap for linux')
class TestNormalPmap(unittest.TestCase):

    def setUp(self):
        self.traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")

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
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")

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

         # test worker
         # need to test this since coverages seems not recognize partial func
        data = worker_byfunc(rank=2, n_cores=8, func=pt.radgyr, traj=traj, args=(), kwd={'mask': '@CA'}, iter_options={})
        assert data[0] == 2, 'rank must be 2'
        assert data[2] == 1, 'n_frames for rank=2 should be 1 (only 10 frames in total)'

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

@unittest.skipUnless(sys.platform.startswith('linux'), 'pmap for linux')
class TestParallelMapForTrajectoryIteratorWithTransformation(unittest.TestCase):
    def test_trajectoryiterator_with_transformation(self):
        traj_on_disk = pt.iterload("data/tz2.nc", "data/tz2.parm7")
        traj_on_mem  = pt.load("data/tz2.nc", "data/tz2.parm7")

        traj_on_disk.superpose(mask='@CA', ref=3)
        traj_on_mem.superpose(mask='@CA', ref=3)

        rmsd0_dict = pt.pmap(pt.rmsd_nofit, traj_on_disk, mask='@CB', n_cores=2, ref=0)
        rmsd1 = pt.rmsd_nofit(traj_on_mem, mask='@CB', ref=0)

        aa_eq(pt.tools.dict_to_ndarray(rmsd0_dict),
              rmsd1)

@unittest.skipUnless(sys.platform.startswith('linux'), 'pmap for linux')
class TestParallelMapForMatrix(unittest.TestCase):

    def test_matrix_module(self):
        traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")

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


@unittest.skipUnless(sys.platform.startswith('linux'), 'pmap for linux')
class TestCpptrajCommandStyle(unittest.TestCase):

    def test_c_command_style(self):
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
            # another way to get reference
            data2 = pt.pmap(['rms @CA reference'],
                            traj,
                            ref=traj[3],
                            n_cores=n_cores)
            # use int for ref
            data3 = pt.pmap(pt.rmsd,
                            traj,
                            ref=3,
                            mask='@CA',
                            n_cores=n_cores)
            # use int for ref: use cpptraj's commmand style
            data4 = pt.pmap(['rms @CA reference'],
                            traj,
                            ref=3,
                            n_cores=n_cores)
            arr = pt.tools.dict_to_ndarray(data)[0]
            arr2 = pt.tools.dict_to_ndarray(data2)[0]
            arr3 = pt.tools.dict_to_ndarray(data3)[0]
            arr4 = pt.tools.dict_to_ndarray(data4)[0]
            aa_eq(arr, pt.rmsd(traj, ref=3, mask='@CA'))
            aa_eq(arr2, pt.rmsd(traj, ref=3, mask='@CA'))
            aa_eq(arr3, pt.rmsd(traj, ref=3, mask='@CA'))
            aa_eq(arr4, pt.rmsd(traj, ref=3, mask='@CA'))

            # use 4-th and 5-th Frame for reference
            data = pt.pmap(
                ['rms @CA refindex 0', 'rms @CB refindex 1'],
                traj,
                ref=[traj[3], traj[4]],
                n_cores=n_cores)
            arr = pt.tools.dict_to_ndarray(data)[0]
            aa_eq(arr, pt.rmsd(traj, '@CA', 3))

            arr1 = pt.tools.dict_to_ndarray(data)[1]
            aa_eq(arr1, pt.rmsd(traj, ref=4, mask='@CB'))

            # no ref
            data = pt.pmap(['radgyr', ], traj, n_cores=n_cores)
            arr = pt.tools.dict_to_ndarray(data)[0]
            aa_eq(arr, pt.radgyr(traj))


@unittest.skipUnless(sys.platform.startswith('linux'), 'pmap for linux')
class TestParallelMapForAverageStructure(unittest.TestCase):

    def test_pmap_average_structure(self):
        traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")
        saved_frame = pt.mean_structure(traj, '@CA')
        saved_xyz = saved_frame.xyz

        for n_cores in [2, 3, 4]:
            frame = pt.pmap(pt.mean_structure, traj, '@CA', n_cores=n_cores)
            aa_eq(frame.xyz, saved_xyz)


@unittest.skipUnless(sys.platform.startswith('linux'), 'pmap for linux')
class TestLoadBathPmap(unittest.TestCase):

    def test_load_batch(self):
        '''just test ValueError
        '''
        self.assertRaises(
            ValueError,
            lambda: _load_batch_pmap(n_cores=4, lines=['autoimage'], traj=None, dtype='dict', root=0, mode='xyz', ref=None))


@unittest.skipUnless(sys.platform.startswith('linux'), 'pmap for linux')
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
                parallel_out_c_style = pt.pmap(
                    ['radgyr @CA nomax'],
                    traj,
                    frame_indices=frame_indices)
                aa_eq(serial_out, pt.tools.dict_to_ndarray(parallel_out))
                aa_eq(serial_out,
                      pt.tools.dict_to_ndarray(parallel_out_c_style))


@unittest.skipUnless(sys.platform.startswith('linux'), 'pmap for linux')
class TestCheckValidCommand(unittest.TestCase):

    def test_check_valid_command(self):
        from pytraj.parallel.base import check_valid_command
        assert check_valid_command(['rms', ]) == (['rms refindex 0 '], True)
        assert check_valid_command(['distrmsd', ]) == (['distrmsd refindex 0 '], True)
        assert check_valid_command(['nativecontacts', ]) == (['nativecontacts refindex 0 '], True)
        assert check_valid_command(['nastruct', ]) == (['nastruct refindex 0 '], True)
        assert check_valid_command(['symmetricrmsd', ]) == (['symmetricrmsd refindex 0 '], True)
        traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")

        aa_eq(pt.tools.dict_to_ndarray(
            pt.pmap(['rmsd'], traj, ref=traj[3], n_cores=3)),
            pt.rmsd(traj, ref=traj[3]))

        # provide refindex
        aa_eq(pt.tools.dict_to_ndarray(
            pt.pmap(['rmsd refindex 0'], traj, ref=traj[3], n_cores=3)),
            pt.rmsd(traj, ref=traj[3]))

        aa_eq(pt.tools.dict_to_ndarray(
            pt.pmap(['rmsd refindex 0'], traj, ref=[traj[3], traj[0]], n_cores=3)),
            pt.rmsd(traj, ref=traj[3]))

        # if user does not provide reference, need to give it to them
        aa_eq(pt.tools.dict_to_ndarray(
            pt.pmap(['rmsd'], traj, n_cores=3)),
            pt.rmsd(traj, ref=traj[0]))

        # does not support matrix
        self.assertRaises(ValueError, lambda: pt.pmap(['matrix'], traj, n_cores=2))

        # do not accept any c analysis command
        for word in c_commands.analysis_commands:
            self.assertRaises(ValueError, lambda: pt.pmap(word, traj, n_cores=2))


@unittest.skipUnless(sys.platform.startswith('linux'), 'pmap for linux')
class TestVolmap(unittest.TestCase):

    def test_volmap(self):
        traj = pt.iterload("data/tz2.ortho.nc", "data/tz2.ortho.parm7")

        # raise if does not provide size
        self.assertRaises(AssertionError, lambda: pt.pmap(pt.volmap, traj, mask=':WAT@O',
                                                          grid_spacing=(0.5, 0.5, 0.5),
                                                          n_cores=2))

        mask = ':WAT@O'
        grid_spacing = (0.5, 0.5, 0.5)

        for n_cores in [1, 2, 3]:
            for size in [(20, 20, 20), (20, 40, 60)]:
                serial_out = pt.volmap(traj, mask=mask, grid_spacing=grid_spacing, size=size)
                parallel_out = pt.pmap(pt.volmap, traj, mask=mask, grid_spacing=grid_spacing,
                                       size=size, n_cores=n_cores)
                self.assertEqual(serial_out.shape, tuple(2 * x for x in size))
                aa_eq(serial_out, parallel_out)


@unittest.skipUnless(sys.platform.startswith('linux'), 'pmap for linux')
class TestWorker(unittest.TestCase):

    def testworker_by_actlist(self):
        # just want to exercise all codes
        traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")
        for ref in [None, traj[0], [traj[0], traj[1]]]:
            data = worker_by_actlist(rank=3, n_cores=8, traj=traj, lines=['radgyr @CA', 'vector :3 :7'],
                                     ref=ref, kwd=dict())


def change_10_atoms(traj):
    for frame in traj:
        frame.xyz[:10] += 1.
        yield frame


@unittest.skipUnless(sys.platform.startswith('linux'), 'pmap for linux')
class TestInserNewFunction(unittest.TestCase):

    def test_insert_new_function(self):
        traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")

        # create mutable Trajectory
        t0 = traj[:]
        for frame in t0:
            frame.xyz[:10] += 1.

        data_parallel = pt.pmap(pt.radgyr, traj, n_cores=2, apply=change_10_atoms)
        data_serial = pt.radgyr(t0)
        aa_eq(data_parallel['RoG_00000'], data_serial)


if __name__ == "__main__":
    unittest.main()
