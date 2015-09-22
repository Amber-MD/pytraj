from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.testing import cpptraj_test_dir, duplicate_traj
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        # all atoms
        fa = mdio.load(traj(0, 5), traj.top)
        fa2 = traj[:6]
        aa_eq(fa.xyz, fa2.xyz)

        # strip atoms
        new_top = traj.top.strip_atoms("!@CA", copy=True)
        fa = mdio.load(traj(0, 5, mask='@CA'), top=new_top)
        fa2 = traj[:6, '@CA']
        aa_eq(fa.xyz, fa2.xyz)


if __name__ == "__main__":
    unittest.main()
