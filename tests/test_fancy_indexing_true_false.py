from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):

    def test_0(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        # mutable traj
        traj2 = traj.to_mutable_trajectory()
        # NotImplementedError: mutable traj
        self.assertRaises(NotImplementedError, lambda: traj2[[True, False]])
        # NotImplementedError : trajiter
        self.assertRaises(NotImplementedError, lambda: traj[[True, False]])

if __name__ == "__main__":
    unittest.main()
