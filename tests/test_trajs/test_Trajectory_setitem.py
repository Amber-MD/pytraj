from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
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

        # a set of atoms, two trajectories
        fa0 = traj[:]
        fa1 = traj[:]
        fa1 += 1.
        aa_eq(fa0.xyz + 1., fa1.xyz)
        atm = traj.top("!@H=")
        # try assigning from AtomMask
        fa0[atm] = fa1[atm]
        aa_eq(fa0[atm].xyz, fa1[atm].xyz)
        aa_eq(fa0['@H='].xyz, traj['@H='].xyz)

        # all atoms
        xyz = traj.xyz + 2.
        fa["*"] = xyz
        aa_eq(fa.xyz, xyz)

        # all atoms for a set of frames
        xyz = traj.xyz[:3] + 3.
        fa[:3]["*"] = xyz
        aa_eq(fa.xyz[:3], xyz)

        # automatically cast
        fa0 = fa.copy()
        xyz = fa.xyz + 1.
        fa0[0] = xyz[0]  # fa[0] return a Frame
        aa_eq(fa0[0].xyz, xyz[0])
        # try to assign a Frame
        #print(fa0, fa)
        fa0[0] = fa[0]
        aa_eq(fa0[0].xyz, fa[0].xyz)

        def shape_mismatch():
            fa[0] = xyz

        self.assertRaises(ValueError, lambda: shape_mismatch())

        def shape_mismatch2():
            fa[0] = Frame()

        self.assertRaises(ValueError, lambda: shape_mismatch2())

        # assign to None
        def None_value():
            fa[0] = None

        self.assertRaises(ValueError, lambda: None_value())


if __name__ == "__main__":
    unittest.main()
