from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca
from pytraj.utils import Timer

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        fa = traj[:]
        print (fa)
        atm = traj.top("@CA")
        indices = atm.indices
        xyz = traj['@CA'].xyz.copy() + 1.
        X = traj.xyz.copy()

        assert X[:, indices].shape == xyz.shape
        assert fa['@CA'].shape == xyz.shape

        # assignment speed
        # 10 frames
        @Timer()
        def speed_test_np():
            X[:, indices] = xyz

        @Timer()
        def speed_test_pytraj():
            fa['@CA'] = xyz

        speed_test_np() # 3 times faster
        speed_test_pytraj()
        # make sure ther are equal
        aa_eq(fa['@CA'].xyz, X[:, indices])

        # assignment speed
        # 50K frames
        FA = fa.copy()
        for _ in range(4):
            FA += FA + FA + FA
        FA += FA[:10000].copy()
        fa = FA[:50000]
        print (fa)
        xyz = fa['@CA'].xyz.copy() + 1.
        X = fa.xyz.copy()
        @Timer()
        def speed_test_np():
            X[:, indices] = xyz

        @Timer()
        def speed_test_pytraj():
            fa['@CA'] = xyz

        speed_test_np() # 
        speed_test_pytraj() # similiar speed
        # make sure ther are equal
        aa_eq(fa['@CA'].xyz, X[:, indices])

        # assignment speed
        # 50K frames, different mask
        mask = "@H=,CB"
        indices = traj.top(mask).indices
        del xyz
        xyz = fa[mask].xyz.copy() + 1.
        X = fa.xyz.copy()
        @Timer()
        def speed_test_np():
            X[:, indices] = xyz

        @Timer()
        def speed_test_pytraj():
            fa[mask] = xyz

        speed_test_np() # 
        speed_test_pytraj() # similiar speed
        # make sure ther are equal
        aa_eq(fa[mask].xyz, X[:, indices])

        # assignment speed
        # 50K frames, different mask
        mask = "!@H=,CB"
        indices = traj.top(mask).indices
        del xyz
        xyz = fa[mask].xyz.copy() + 1.
        X = fa.xyz.copy()
        @Timer()
        def speed_test_np():
            X[:, indices] = xyz

        @Timer()
        def speed_test_pytraj():
            fa[mask] = xyz

        speed_test_np() # 
        speed_test_pytraj() # similiar speed
        # make sure ther are equal
        aa_eq(fa[mask].xyz, X[:, indices])


if __name__ == "__main__":
    unittest.main()
