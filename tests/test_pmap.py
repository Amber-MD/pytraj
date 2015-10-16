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


class TestParallelMapForMatrix(unittest.TestCase):
    def test_matrices(self):
        traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")

        # not support [covar, distcovar, mwcovar]
        for func in [matrix.dist, matrix.idea]:
            x = pt.pmap(4, func, traj, '@CA')
            y = np.sum((val[1] * val[2] for val in x), axis=1)
            aa_eq(y/traj.n_frames, func(traj, '@CA'))


if __name__ == "__main__":
    unittest.main()
