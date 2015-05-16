from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq, eq_coords
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir, duplicate_traj
import pytraj.common_actions as pyca
from pytraj.utils import Timer

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        fa = traj[:]

        @Timer()
        def test_time_aa_eq():
            aa_eq(fa.xyz, traj.xyz)

        @Timer()
        def test_time_eq_coords():
            eq_coords(fa, traj)

        print ("test_time_aa_eq")
        test_time_aa_eq()
        print ("test_time_eq_coords")
        test_time_eq_coords()

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
        # 10000 frames
        FA = duplicate_traj(traj, 1000)
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
        #aa_eq(fa['@CA'].xyz, X[:, indices])

        # assignment speed
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

        # assignment speed
        xyz = fa.xyz + 10.
        fa2 = fa.copy()

        @Timer()
        def update_xyz_version():
            fa.update_xyz(xyz)

        @Timer()
        def fancy_indexing_version():
            fa2['*'] = xyz

        print ("update_xyz_version")
        update_xyz_version()
        print ("fancy_indexing_version")
        fancy_indexing_version()
        aa_eq(fa.xyz ,fa2.xyz) # ok


if __name__ == "__main__":
    unittest.main()
