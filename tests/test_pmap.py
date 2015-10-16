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

    def test_regular1D(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        func_list = [pt.radgyr, pt.molsurf, pt.rmsd]
        ref = traj[-3]

        for n_cores in [2, 3, 4]:
            for func in func_list:
                if func in [pt.rmsd, ]:
                    pout = gather(pt.pmap(n_cores, func, traj, ref=ref))
                    serial_out = flatten(func(traj, ref=ref))
                else:
                    pout = gather(pt.pmap(n_cores, func, traj))
                    serial_out = flatten(func(traj))
                aa_eq(pout, serial_out)

        # search_hbonds
        a = pt.pmap(4, pt.search_hbonds, traj, dtype='dataset')
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
            pt.pmap(2, pt.bfactors, traj[:])

        self.assertRaises(ValueError, lambda: need_to_raise_2())

    def test_different_references(self):
        traj = self.traj
        func = pt.rmsd
        for i in range(0, 8, 2):
            ref = self.traj[i]
            for n_cores in [2, 3, 4, 5]:
                pout = gather(pt.pmap(n_cores, func, traj, ref=ref))
                serial_out = flatten(func(traj, ref=ref))
                aa_eq(pout, serial_out)

class TestParallelMapForMatrix(unittest.TestCase):
    def test_matrices(self):
        traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")

        # not support [covar, distcovar, mwcovar]
        for n_cores in [2, 3, 4, 5]:
            for func in [matrix.dist, matrix.idea]:
                x = pt.pmap(n_cores, func, traj, '@CA')
                y = np.sum((val[1] * val[2] for val in x))
                aa_eq(y/traj.n_frames, func(traj, '@CA'))

class TestCpptrajCommandStyle(unittest.TestCase):
    def test_cpptraj_command_style(self):
        traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")

        angle_ = pt.angle(traj, ':3 :4 :5')
        distance_ = pt.distance(traj, '@10 @20')

        data = pt.pmap(4, ['angle :3 :4 :5', 'distance @10 @20'], traj)
        aa_eq(angle_, data['Ang_00002'])
        aa_eq(distance_, data['Dis_00003'])

class TestParallelMapForAverageStructure(unittest.TestCase):
    def test_pmap_average_structure(self):
        traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")
        saved_frame = pt.mean_structure(traj, '@CA')
        saved_xyz = saved_frame.xyz

        for n_cores in [2, 3, 4, 5]:
            data = pt.pmap(n_cores, pt.mean_structure, traj, '@CA')
            xyz = np.sum((x[2] * np.array(x[1].xyz) for x in data)) / traj.n_frames
            aa_eq(xyz, saved_xyz)


if __name__ == "__main__":
    unittest.main()
