from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class Test(unittest.TestCase):

    def test_0(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        # rmsd
        flist = [pt.rmsd,
                 pt.radgyr,
                 pt.molsurf,
                 pt.calc_atomicfluct,
                 pt.bfactors,
                 # pt.calc_rmsd_with_rotation_matrices,
                 pt.calc_pairwise_rmsd, ]
        for func in flist:
            aa_eq(
                pt.tools.flatten(func(traj,
                                      mask=range(7))),
                pt.tools.flatten(func(traj,
                                      mask="@1,2,3,4,5,6,7")))

            aa_eq(
                pt.tools.flatten(func(traj,
                                      mask=range(0, 7, 2))),
                pt.tools.flatten(func(traj,
                                      mask="@1,3,5,7")))


if __name__ == "__main__":
    unittest.main()
