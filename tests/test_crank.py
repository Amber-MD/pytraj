from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):

    def test_0(self):
        # TODO: assert
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        d0 = pt.distance(traj, [0, 50])[0]
        d1 = pt.distance(traj, [10, 100])[0]
        d = pt.crank(d0, d1, mode='distance')

if __name__ == "__main__":
    unittest.main()
