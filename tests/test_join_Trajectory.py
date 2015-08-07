from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        fa0 = traj.copy()
        fa1 = traj.copy()

        # +=
        SIZE = 10
        fa0.join(fa1, copy=False)
        assert fa0.size == traj.size * 2 == 20
        aa_eq(fa0[:SIZE].xyz, fa0[SIZE:].xyz)
        aa_eq(fa0[:SIZE].xyz, traj.xyz)

        # memview
        print(fa0)
        print(fa0[SIZE])
        fa1[0, 0, 0] = 1000.
        assert fa1[0, 0, 0] == 1000.
        assert fa0[SIZE, 0, 0] == 1000.


if __name__ == "__main__":
    unittest.main()
