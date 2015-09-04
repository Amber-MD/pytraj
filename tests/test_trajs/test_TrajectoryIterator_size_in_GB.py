from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        print(traj._estimated_GB)
        traj._size_limit_in_MB = traj._estimated_GB - 0.1
        self.assertRaises(MemoryError, lambda: traj.xyz)

        traj._force_load = True
        traj.xyz


if __name__ == "__main__":
    unittest.main()
