from __future__ import print_function
import unittest
import numpy as np
import pytraj as pt
from pytraj.utils import eq, aa_eq
import pytraj.common_actions as pyca
from pytraj.tools import flatten
from pytraj import matrix
from pytraj.compat import set


def gather(pmap_out):
    pmap_out = sorted(pmap_out, key=lambda x: x[0])
    return flatten([x[1] for x in pmap_out])


class TestNormal(unittest.TestCase):
    def setUp(self):
        self.traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

    def test_raise_if_func_not_callable(self):
        self.assertRaises(ValueError, lambda: pt.pmap('test', self.traj))

    def test_regular1D(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        func_list = [pt.radgyr, pt.molsurf, pt.rmsd]
        ref = traj[-3]

        for n_cores in [2, 3, 4]:
            for func in func_list:
                if func in [pt.rmsd, ]:
                    pout = gather(pt.pmap(n_cores=n_cores, func=func, traj=traj, ref=ref))
                    serial_out = flatten(func(traj, ref=ref))
                else:
                    pout = gather(pt.pmap(n_cores=n_cores, func=func, traj=traj))
                    serial_out = flatten(func(traj))
                aa_eq(pout, serial_out)

        # search_hbonds
        a = pt.pmap(pt.search_hbonds, traj, dtype='dataset', n_cores=4)
        pout = pt.tools.flatten([x[1]['total_solute_hbonds'] for x in a])
        serial_out = pt.search_hbonds(traj, dtype='dataset')['total_solute_hbonds']
        aa_eq(pout, serial_out)

        keys = pt.tools.flatten([x[1].keys() for x in a])

        # raise if a given method does not support pmap
        def need_to_raise(traj=traj):
            pt.pmap(2, pt.bfactors, traj)

        self.assertRaises(ValueError, lambda: need_to_raise())

        # raise if a traj is not TrajectoryIterator
        def need_to_raise_2(traj=traj):
            pt.pmap(pt.bfactors, traj[:], n_cores=2)

        self.assertRaises(ValueError, lambda: need_to_raise_2())

    def test_different_references(self):
        traj = self.traj
        func = pt.rmsd
        for i in range(0, 8, 2):
            ref = self.traj[i]
            for n_cores in [2, 3, 4, 5]:
                pout = gather(pt.pmap(n_cores=n_cores, func=func, traj=traj, ref=ref))
                serial_out = flatten(func(traj, ref=ref))
                aa_eq(pout, serial_out)

    def test_iter_options(self):
        traj = pt.iterload("data/tz2.ortho.nc", "data/tz2.ortho.parm7")
        t0 = traj[:].autoimage().rmsfit(ref=0)
        saved_avg = pt.mean_structure(t0)

        # perform autoimage, then rms fit to 1st frame, then compute mean structure
        iter_options = {'autoimage': True, 'rmsfit': 0}
        for n_cores in [2, 3]:
            avg = pt.pmap(pt.mean_structure, traj, iter_options=iter_options,
                    n_cores=n_cores)
            aa_eq(saved_avg.xyz, avg.xyz)


class TestParallelMapForMatrix(unittest.TestCase):
    def test_matrix_module(self):
        traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")

        # not support [covar, distcovar, mwcovar]
        for n_cores in [2, 3, 4, 5]:
            for func in [matrix.dist, matrix.idea]:
                x = pt.pmap(func, traj, '@CA', n_cores=n_cores)
                aa_eq(x, func(traj, '@CA'))

    def test_ired_vector_and_matrix_pmap(self):
        traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")
        h = traj.top.select('@H')
        n = h - 1
        nh = list(zip(n ,h))

        exptected_vecs, exptected_mat = pt.ired_vector_and_matrix(traj, nh)
        for n_cores in [2, 4, 6]:
            vecs, mat = pt.pmap(pt.ired_vector_and_matrix, traj, nh, n_cores=n_cores)
            aa_eq(exptected_vecs, vecs, decimal=7)
            aa_eq(exptected_mat, mat, decimal=7)

    def test_rotation_matrix_in_rmsd_calculation(self):
        traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")
        saved_mat = pt.rotation_matrix(traj, ref=traj[3], mask='@CA')
        saved_rmsd = pt.rmsd(traj, ref=traj[3], mask='@CA')

        for n_cores in [2, 3, 4, 5]:
            mat = pt.pmap(pt.rotation_matrix, traj, ref=traj[3], mask='@CA')
            mat2, rmsd_  = pt.pmap(pt.rotation_matrix, traj, ref=traj[3], mask='@CA',
                    with_rmsd=True)
            aa_eq(saved_mat, mat)
            aa_eq(saved_mat, mat2)
            aa_eq(saved_rmsd, rmsd_)

class TestCpptrajCommandStyle(unittest.TestCase):
    def test_cpptraj_command_style(self):
        traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")

        angle_ = pt.angle(traj, ':3 :4 :5')
        distance_ = pt.distance(traj, '@10 @20')

        data = pt.pmap(['angle :3 :4 :5', 'distance @10 @20'], traj, n_cores=2)
        aa_eq(angle_, data['Ang_00002'])
        aa_eq(distance_, data['Dis_00003'])

class TestParallelMapForAverageStructure(unittest.TestCase):
    def test_pmap_average_structure(self):
        traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")
        saved_frame = pt.mean_structure(traj, '@CA')
        saved_xyz = saved_frame.xyz

        for n_cores in [2, 3, 4, 5]:
            frame = pt.pmap(pt.mean_structure, traj, '@CA', n_cores=n_cores)
            aa_eq(frame.xyz, saved_xyz)


if __name__ == "__main__":
    unittest.main()
