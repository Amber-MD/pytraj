from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
import pytraj.common_actions as pyca
from pytraj.tools import flatten


def gather(pmap_out):
    pmap_out = sorted(pmap_out, key=lambda x: x[0])
    return flatten([x[1] for x in pmap_out])


class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        # n_cores = 3
        # radgyr
        # TODO: hang forever with pt.rmsd
        #func_list = [pt.radgyr, pt.molsurf, pt.rmsd]
        func_list = [pt.radgyr, pt.molsurf]
        ref = traj[-3]

        for n_cores in [2, 3, 4]:
            for func in func_list:
                if func in [pt.rmsd, ]:
                    print(func)
                    pout = gather(pt.pmap(n_cores, func, traj, ref=ref))
                    serial_out = flatten(func(traj, ref=ref))
                else:
                    pout = gather(pt.pmap(n_cores, func, traj))
                    serial_out = flatten(func(traj))
                aa_eq(pout, serial_out)
        # search_hbonds
        a = pt.pmap(4, pt.search_hbonds, traj)
        pout = pt.tools.flatten([x[1]['total_solute_hbonds'] for x in a])
        serial_out = pt.search_hbonds(traj)['total_solute_hbonds']
        #print(pout, serial_out)
        aa_eq(pout, serial_out)

        keys = pt.tools.flatten([x[1].keys() for x in a])
        from pytraj.compat import set

        #print(set(keys))

        # raise if a given method does not support pmap
        def need_to_raise(traj=traj):
            pt.pmap(2, pt.bfactors, traj)

        self.assertRaises(ValueError, lambda: need_to_raise())

        # raise if a traj is not TrajectoryIterator
        def need_to_raise_2(traj=traj):
            pt.pmap(2, pt.bfactors, traj[:])

        self.assertRaises(ValueError, lambda: need_to_raise_2())

        #need_to_raise()


if __name__ == "__main__":
    unittest.main()
