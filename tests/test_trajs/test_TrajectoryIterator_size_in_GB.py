from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        print(traj._estimated_GB)
        traj._size_limit_in_GB = traj._estimated_GB - 0.1
        self.assertRaises(MemoryError, lambda: traj.xyz)

        traj._force_load = True
        traj.xyz


if __name__ == "__main__":
    unittest.main()
