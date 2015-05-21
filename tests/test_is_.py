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

        # test is_
        assert fa[0].is_(fa[0])
        assert not fa[0].is_(fa[1])

        fa_sliced = fa[3:7]
        assert fa_sliced[0].is_(fa[3])

        # test __contains__
        f0 = fa[0]
        assert f0 in fa

        f0_from_iter = traj[0]
        assert not f0_from_iter in fa

        # make a new copy
        f1 = f0.copy()
        aa_eq(f1.xyz, f0.xyz)
        assert not f1 in fa

        # append
        fa.append(f1, copy=False)
        assert f1 in fa

        # append
        f2 = fa[5] + 1.
        fa.append(f2, copy=True)
        assert not f2 in fa

if __name__ == "__main__":
    unittest.main()
