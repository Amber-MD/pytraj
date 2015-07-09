from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.testing import aa_eq
from pytraj.utils.check_and_assert import is_word_in_class_name
from pytraj.decorators import no_test, test_if_having
from pytraj.compat import zip
from pytraj.utils import Timer
from pytraj.core.TrajinList import TrajinList

fname = "./data/md1_prod.Tc5b.x"
topname = "./data/Tc5b.top"
traj = mdio.iterload(fname, topname)
dslist = DataSetList()
dslist.add_set("coords")
dslist[0].top = traj.top.copy()
dslist[0].load(fname)
dslist.add_set("traj")
dslist[1].top = traj.top.copy()
dslist[1].load(fname)
farray = traj[:]
tlist = TrajinList()
tlist.top = traj.top.copy()
tlist.add_traj(fname)


@Timer()
def test_load(my_traj, ref_traj=None, n_frames=None):
    print()
    if ref_traj is None:
        ref_traj = traj

    if n_frames is None:
        n_frames = ref_traj.n_frames

    fa = Trajectory()
    fa.top = traj.top.copy()
    fa.load(my_traj)

    assert fa.size == n_frames

    for f0, f1 in zip(fa, ref_traj):
        assert_almost_equal(f0.coords, f1.coords)


class Test(unittest.TestCase):

    @test_if_having("numpy")
    def test_0(self):
        import numpy as np
        print("load from TrajectoryIterator")
        test_load(traj)
        print("load from Trajectory")
        test_load(farray)
        print(" load from frame_iter, TrajectoryIterator")
        test_load(traj())
        print(" load from frame_iter, Trajectory")
        test_load(farray())
        print(" load from file")
        test_load(fname)
        print(" load from a list of files")
        test_load([fname for _ in range(10)], n_frames=traj.n_frames * 10)
        print(" load from TrajectoryIterator")
        test_load(traj)
        print(" load from DataSet Coords")
        test_load(dslist[0])
        #print (" load from DataSet Traj, sometimes got segmentation fault")
        # test_load(dslist[1])
        print(" load from ndarray")
        test_load(traj.xyz)
        print(" load from ndarray, flatten")
        test_load(traj.xyz.flatten())
        print(" load from ndarray, flatten")
        test_load(traj.xyz.flatten().tolist())
        print(" load from ndarray, flatten")
        test_load(traj.xyz.tolist())

    @test_if_having("mdtraj")
    def test_1(self):
        import mdtraj as md
        import numpy as np

        traj = mdio.iterload(fname, topname)

        m_top = md.load_prmtop(topname)
        m_traj = md.load_mdcrd(fname, m_top)
        print(m_traj.xyz[0, 0, 0])
        print(traj[0, 0, 0])
        indices = np.empty((20, 2), dtype=np.int64)
        indices[:, 0] = traj.top("@CA").selected_indices()[:20]
        indices[:, 1] = traj.top("@H*").selected_indices()[:20]

        _fa = Trajectory()
        _fa.top = traj.top.copy()
        _fa.load(m_traj)  # not auto-cast from `nm` to `angstrom`
        aa_eq(_fa.xyz * 10., traj.xyz)

        a_mdtraj = md.compute_distances(m_traj, indices)
        a_fa_from_mdtraj = _fa.calc_distance(indices)
        import pytraj.common_actions as pyca
        a_traj = pyca.calc_distance(traj, indices)
        aa_eq(a_mdtraj,  a_fa_from_mdtraj)
        aa_eq(a_mdtraj * 10,  a_traj)
        aa_eq(a_fa_from_mdtraj * 10,  a_traj)

        # load mdtraj with given Topology
        fa2 = Trajectory(m_traj, top=_fa.top)
        aa_eq(fa2.xyz, _fa.xyz)

    @no_test
    def test_2(self):
        # turn off this test since getting 2/3 chances of segmentation fault
        # Don't know why
        print("test loading DataSetList")
        _dslist = dslist
        print(dslist)
        ref_traj = Trajectory()
        ref_traj.top = traj.top.copy()
        ref_traj.load(traj)
        ref_traj.load(traj)
        test_load(_dslist, ref_traj=ref_traj)

    def test_3(self):
        print("test loading frame")
        fa = Trajectory(traj[0], traj.top)
        assert fa.size == 1
        assert_almost_equal(fa[0].coords, traj[0].coords)

if __name__ == "__main__":
    unittest.main()
