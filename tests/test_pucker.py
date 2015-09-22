from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.load_pdb("./data/Test_NAstruct/adh026.3.pdb")
        d = pt.common_actions.pucker(traj, resrange=range(3, 7))


if __name__ == "__main__":
    unittest.main()
