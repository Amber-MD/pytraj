from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca
from pytraj import AtomMask

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        fa = traj[:]

        # single value
        fa[0, 0, 0] = 100.
        assert fa[0, 0, 0] == 100.

        # single atom
        fa[0, 0] = [100., 101, 102]
        aa_eq(fa[0, 0], [100, 101, 102])

        # a set of atoms
        indices = [1, 10, 11]
        atm = AtomMask(indices)
        xyz_sub = fa.xyz[:, indices] + 1.
        fa[atm] = xyz_sub
        aa_eq(fa[atm].xyz, xyz_sub)

        indices = traj.top("@CA")
        atm = AtomMask(indices)
        xyz_sub = fa.xyz[:, list(indices)] + 1.
        fa[atm] = xyz_sub
        aa_eq(fa['@CA'].xyz, xyz_sub)

        # all atoms
        xyz = traj.xyz + 2.
        fa["*"] = xyz
        aa_eq(fa.xyz, xyz)

        # all atoms for a set of frames
        xyz = traj.xyz[:3] + 3.
        fa[:3]["*"] = xyz
        aa_eq(fa.xyz[:3], xyz)

if __name__ == "__main__":
    unittest.main()
