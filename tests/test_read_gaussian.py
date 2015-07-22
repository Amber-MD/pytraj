from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):

    @test_if_having("cclib")
    def test_0(self):
        # TODO : assert
        traj = pt.tools.read_gaussian_output("./data/gaussian/GF2.log",
                                             "./data/gaussian/GF2.pdb")
        print (traj)

if __name__ == "__main__":
    unittest.main()
