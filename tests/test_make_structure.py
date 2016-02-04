from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestMakeStructure(unittest.TestCase):

    def test_makestructure(self):
        # https://github.com/Amber-MD/cpptraj/issues/27
        # load only 1st frame
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")

        #  pply polyproline II dihedral to residues 1-13
        t0 = traj[:1].copy()
        pt.make_structure(t0, 'pp2:1-13')
        pt.write_traj('./output/test0.pdb',
                      t0,
                      options='model',
                      overwrite=True)

        t1 = traj[:1].copy()
        pt.make_structure(t1, "chi1:8:N:CA:CB:CG:35")
        pt.write_traj('./output/test1.pdb',
                      t1,
                      options='model',
                      overwrite=True)


if __name__ == "__main__":
    unittest.main()
