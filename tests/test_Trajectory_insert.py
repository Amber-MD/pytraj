from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir, duplicate_traj
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        fa = traj[:]
        fa2 = traj[:]
        fa3 = traj[:]

        # copy=True
        fa.insert(fa2[0], 7, copy=True)
        assert fa.size == fa2.size + 1
        print(fa[7, 0])
        print(fa2[0, 0])
        assert fa[7].rmsd_nofit(fa2[0]) == 0.
        assert not fa2[0] in fa

        # copy=False
        fa = traj[:]
        fa2 = traj[:]
        fa3 = traj[:]
        fa.insert(fa2[0], 7, copy=False)
        assert fa.size == fa2.size + 1
        print(fa[7, 0])
        print(fa2[0, 0])
        assert fa[7].rmsd_nofit(fa2[0]) == 0.
        assert fa2[0] in fa

        # insert it's self, copy=True
        fa = traj[:]
        fa.insert(fa[0], 5, copy=True)
        assert not fa[0].is_(fa[5])
        fa[0, 0, 0] = 1.5
        assert fa[0, 0, 0] == 1.5
        assert fa[5, 0, 0] != 1.5

        # insert it's self, copy=False
        fa = traj[:]
        fa.insert(fa[0], 5, copy=False)
        assert fa[0].is_(fa[5])
        fa[0, 0, 0] = 1.5
        assert fa[0, 0, 0] == 1.5
        assert fa[5, 0, 0] == 1.5


if __name__ == "__main__":
    unittest.main()
