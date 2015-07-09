from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):

    def test_0(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        # rmsd
        flist = [pt.rmsd, pt.radgyr, pt.molsurf,
                 pt.calc_atomicfluct, pt.bfactors,
                 pt.calc_jcoupling,
                 pt.calc_pairwise_rmsd,
                 pt.calc_rmsd_with_rotation_matrices,
                ]
        for func in flist:
            aa_eq(func(traj, mask=range(7)).flatten(),
                  func(traj, mask="@1,2,3,4,5,6,7").flatten())

            aa_eq(func(traj, mask=range(0, 7, 2)).flatten(),
                  func(traj, mask="@1,3,5,7").flatten())
            print ('%s: OK' % func.__name__)

if __name__ == "__main__":
    unittest.main()
