from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca

def gather(pmap_out):
    from pytraj.tools import flatten
    return flatten([x[1] for x in pmap_out])


class Test(unittest.TestCase):

    def test_0(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        # n_cores = 3
        # radgyr
        func_list = [pt.radgyr, pt.molsurf]
        for n_cores in [2, 3, 4]:
            for func in func_list:
                pout = gather(pt.pmap(n_cores, func, traj))
                print(pout)
                serial_out = func(traj)
                aa_eq(pout, serial_out)

        # reference
        ref = traj[0]
        pout = pt.pmap(4, pt.native_contacts, traj, ref=ref)
        print(pout)

if __name__ == "__main__":
    unittest.main()
