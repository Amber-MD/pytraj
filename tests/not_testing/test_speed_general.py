from __future__ import print_function
import unittest; import pytraj as pt
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

        #print("test_time_aa_eq")
        test_time_aa_eq()
        #print("test_time_eq_coords")
        test_time_eq_coords()

        #print(fa)
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

        speed_test_np()  # 3 times faster
        speed_test_pytraj()
        # make sure ther are equal
        aa_eq(fa['@CA'].xyz, X[:, indices])

        # assignment speed
        # 10000 frames
        FA = duplicate_traj(traj, 1000)
        fa = FA[:50000]
        #print(fa)
        xyz = fa['@CA'].xyz.copy() + 1.
        X = fa.xyz.copy()

        @Timer()
        def speed_test_np():
            X[:, indices] = xyz

        @Timer()
        def speed_test_pytraj():
            fa['@CA'] = xyz

        speed_test_np()
        speed_test_pytraj()  # similiar speed
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

        speed_test_np()
        speed_test_pytraj()  # similiar speed
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

        speed_test_np()
        speed_test_pytraj()  # similiar speed
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

        #print("update_xyz_version")
        update_xyz_version()
        #print("fancy_indexing_version")
        fancy_indexing_version()
        aa_eq(fa.xyz, fa2.xyz)  # ok

        # asarray, xyz
        f = fa[0]
        import numpy as np

        @Timer()
        def call_asarray():
            # to save time, just call
            # f.xyz or f.__array__()
            np.asarray(f)

        @Timer()
        def call__array__():
            f.__array__()

        @Timer()
        def call_xyz():
            f.xyz

        @Timer()
        def call_buffer2d():
            f.buffer2d

        @Timer()
        def call_buffer2d_slice():
            f.buffer2d[:]

        #print("call_asarray")
        call_asarray()  # ~10 times slower
        #print("call__array__")
        call__array__()
        #print("call_xyz")
        call_xyz()
        #print("call_buffer2d")
        call_buffer2d()  # 2-3 times faster
        #print("call_buffer2d_slice")
        call_buffer2d_slice()

    def test_1(self):
        from pytraj.testing import duplicate_traj
        #print("test _fast_update_xyz")
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        indices = traj.top("!@H=").indices
        f = traj[0]
        for i in range(4):
            f.append_xyz(f.xyz)
        xyz = f.xyz + 1.

        @Timer()
        def test_numpy_all_indices():
            xyz[:] = xyz

        @Timer()
        def test_pytraj_all_indices():
            f._fast_copy_from_xyz(xyz)

        xyz[indices] += 1.
        new_xyz = xyz[indices]
        newer_xyz = xyz[indices] + 1.

        @Timer()
        def test_numpy_indices():
            new_xyz[:] = newer_xyz

        @Timer()
        def test_pytraj_indices():
            f._fast_copy_from_xyz(newer_xyz, indices)

        #print("numpy: all xyz")
        test_numpy_all_indices()
        #print("pytraj: all xyz")
        test_pytraj_all_indices()  # similiar speed (need to run several times)
        aa_eq(f.xyz, xyz)

        #print('numpy with indices')
        test_numpy_indices()
        #print('pytraj with indices')
        test_pytraj_indices()  # similiar speed (need to run several times)
        aa_eq(f[indices], new_xyz)


if __name__ == "__main__":
    unittest.main()
