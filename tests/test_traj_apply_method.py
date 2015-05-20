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

        # simple, add 1.0 to all coords
        fa = traj[:]
        xyz = traj.xyz[:]
        fa.apply(lambda x: x + 1.)
        xyz += 1.0
        aa_eq(fa.xyz, xyz)
        aa_eq(traj.xyz, xyz - 1.)

        # simple, add 1.0 to all coords of CA atoms
        fa = traj[:]
        xyz = traj['@CA'].xyz[:]
        fa.apply(lambda x: x + 1., indices_or_mask='@CA')
        xyz += 1.0
        aa_eq(fa['@CA'].xyz, xyz)

        # simple, add 1.0 to atoms 0, 1, 2
        fa = traj[:]
        xyz = traj.xyz[:]
        fa.apply(lambda x: x + 1., indices_or_mask=[0, 1, 2])
        for i in range(xyz.shape[0]):
            xyz[i][[0, 1, 2]] += 1.
        aa_eq(fa.xyz, xyz)

        # simple, add xyz0
        fa = traj[:]
        xyz = traj.xyz[:]
        xyz0 = xyz[0].copy() + 100.

        def func(xyz, xyz0):
            return xyz + xyz0

        fa.apply(func, [xyz0,])
        xyz += xyz0
        aa_eq(fa.xyz, xyz)

        fa.apply(func, xyz0)
        xyz += xyz0
        aa_eq(fa.xyz, xyz)

        xyz1 = xyz0[[0, 1, 2]] # atom 0, 1, 2
        fa.apply(func, args=xyz1, indices_or_mask=[0, 1, 2])
        for i in range(xyz.shape[0]):
            xyz[i][[0, 1, 2]] += xyz1
        aa_eq(fa.xyz, xyz)

if __name__ == "__main__":
    unittest.main()
