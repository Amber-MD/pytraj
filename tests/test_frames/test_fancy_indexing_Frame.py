from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq


class Test(unittest.TestCase):

    def test_0(self):
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        f0 = traj[0]
        xyz = f0.xyz

        # mask, AtomMask
        atm = traj.top("@CA")
        indices = atm.indices
        mask = '@CA'
        f0.top = traj.top
        aa_eq(f0[mask], f0[atm])
        aa_eq(f0[mask], xyz[indices])
        aa_eq(f0[atm, 0], xyz[indices][0])
        aa_eq(f0[atm, 0, 0], xyz[indices][0, 0])
        aa_eq(f0[indices][0, 0], xyz[indices][0, 0])
        aa_eq(f0[atm, 1:10, :], xyz[indices][1:10, :])
        aa_eq(f0[atm, :, 0], xyz[indices][:, 0])
        aa_eq(f0[atm, :, :], xyz[indices][:, :])
        aa_eq(f0[atm, 1, :], xyz[indices][1, :])
        aa_eq(f0[atm, 1:10, :], xyz[indices][1:10, :])
        aa_eq(f0['@CA', 0], xyz[indices][0])
        aa_eq(f0['@CA', 0, 0], xyz[indices][0, 0])
        aa_eq(f0['@CA', :, 0], xyz[indices][:, 0])
        aa_eq(f0['@CA', :, :], xyz[indices][:, :])
        aa_eq(f0['@CA', 1, :], xyz[indices][1, :])

        # if not string_types or AtomMask: behave like normal numpy array
        aa_eq(f0[0], xyz[0])
        aa_eq(f0[:], xyz[:])
        aa_eq(f0[:10], xyz[:10])
        # x, y, z-coords
        aa_eq(f0[:, 0], xyz[:, 0])
        aa_eq(f0[:, 1], xyz[:, 1])
        aa_eq(f0[:, 2], xyz[:, 2])
        aa_eq(f0[20:100, 0], xyz[20:100, 0])
        aa_eq(f0[20:100, 1], xyz[20:100, 1])
        aa_eq(f0[20:100, 2], xyz[20:100, 2])


if __name__ == "__main__":
    unittest.main()
