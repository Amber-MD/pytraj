from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):

    def test_0(self):
        # TODO: assert
        # https://github.com/Amber-MD/cpptraj/issues/27
        # load only 1st frame
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        print(traj)

        #  pply polyproline II dihedral to residues 1-13
        t0 = traj[:1].copy()
        pt.make_structure(t0, 'pp2:1-13')
        pt.write_traj('./output/test0.pdb', t0, mode='model', overwrite=True)

        t1 = traj[:1].copy()
        pt.make_structure(t1,"chi1:8:N:CA:CB:CG:35")
        pt.write_traj('./output/test1.pdb', t1, mode='model', overwrite=True)

if __name__ == "__main__":
    unittest.main()
