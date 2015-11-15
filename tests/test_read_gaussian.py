from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
import pytraj.common_actions as pyca


class Test(unittest.TestCase):

    def test_0(self):
        # TODO : assert
        traj = pt.tools.read_gaussian_output("./data/gaussian/GF2.log",
                                             "./data/gaussian/GF2.pdb")


if __name__ == "__main__":
    unittest.main()
