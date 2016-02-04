from __future__ import print_function
import unittest
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq


class Test(unittest.TestCase):

    def test_0(self):
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")

        # +=
        f0 = traj[0].copy()
        xyz = f0.xyz.copy()
        f0 += 1.
        aa_eq(f0.xyz, xyz + 1.)

        # -=
        f0 = traj[0].copy()
        xyz = f0.xyz.copy()
        f0 -= 1.
        aa_eq(f0.xyz, xyz - 1.)

        # *=
        f0 = traj[0].copy()
        xyz = f0.xyz.copy()
        f0 *= 2.
        aa_eq(f0.xyz, xyz * 2)

        # /= single value
        f0 = traj[0].copy()
        xyz = f0.xyz.copy()
        f0 /= 2.
        aa_eq(f0.xyz, xyz / 2)

        # /= other Frame
        f0 = traj[0].copy()
        xyz = f0.xyz.copy()
        f1 = traj[1].copy()
        xyz1 = f1.xyz.copy()
        f0 /= f1
        aa_eq(f0.xyz, xyz / xyz1)

        # +
        f0 = traj[0].copy()
        xyz = f0.xyz.copy()
        f1 = f0 + 1.
        aa_eq(f1.xyz, xyz + 1.)

        # -
        f0 = traj[0].copy()
        xyz = f0.xyz.copy()
        f1 = f0 - 1.
        aa_eq(f1.xyz, xyz - 1.)

        # *
        f0 = traj[0].copy()
        xyz = f0.xyz.copy()
        f1 = f0 * 10.
        aa_eq(f1.xyz, xyz * 10.)

        # /
        f0 = traj[0].copy()
        xyz = f0.xyz.copy()
        f1 = f0 / 10.
        aa_eq(f1.xyz, xyz / 10.)


if __name__ == "__main__":
    unittest.main()
