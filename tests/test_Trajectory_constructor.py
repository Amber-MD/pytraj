from __future__ import print_function
import unittest; import pytraj as pt
import pytraj as pt
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.compat import izip
from pytraj.utils import Timer
from pytraj.decorators import test_if_having, no_test
from itertools import chain
from pytraj import Trajectory
from pytraj.testing import aa_eq

TRAJ = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")


class Test(unittest.TestCase):
    @Timer()
    def test_0(self):
        #print("from TrajectoryIterator")

        with Timer() as t:
            traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        #print(t.time_gap())

        with Timer() as t:
            for _f0 in traj:
                pass
        #print(t.time_gap())

        with Timer() as t:
            farray = Trajectory(traj, traj.top)
        #print(t.time_gap())

        for f0, f1 in izip(farray, traj):
            assert_almost_equal(f0.coords, f1.coords)

    @Timer()
    def test_1(self):
        #print("from Trajectory")
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        _farray = traj[:]
        farray = Trajectory(_farray, traj.top)

        for f0, f1 in izip(farray, _farray):
            assert_almost_equal(f0.coords, f1.coords)

    @Timer()
    @test_if_having("numpy")
    def test_2(self):
        #print("from xyz array")
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        farray_0 = Trajectory(traj.xyz, traj.top)

        for f0, f1 in izip(farray_0, traj):
            assert_almost_equal(f0.coords, f1.coords)

        farray_1 = Trajectory(xyz=traj.xyz, top=traj.top)
        for f0, f1 in izip(farray_1, traj):
            assert_almost_equal(f0.coords, f1.coords)

    @Timer()
    def test_3(self):
        #print("from dataset Coords")
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist = DataSetList()
        dslist.add_set("coords", "", "")
        dslist[0].top = traj.top.copy()
        dslist[0].load(traj)

        farray_0 = Trajectory(dslist[0], traj.top)
        for f0, f1 in izip(farray_0, traj):
            assert_almost_equal(f0.coords, f1.coords)

    @Timer()
    @test_if_having("numpy")
    def test_5(self):
        #print("from list")
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        farray_0 = Trajectory(traj[:, :, :], traj.top)
        for f0, f1 in izip(farray_0, traj):
            assert_almost_equal(f0.coords, f1.coords)

    @Timer()
    def test_6(self):
        #print("from list 0")
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        # unpack 3D list
        xyz_list = traj.tolist()
        #print(xyz_list.__len__())
        assert len(xyz_list) == traj.n_frames
        farray_0 = Trajectory(xyz_list, traj.top)
        assert farray_0.size == traj.size
        for f0, f1 in izip(farray_0, traj):
            assert_almost_equal(f0.coords, f1.coords)

    @Timer()
    def test_7(self):
        #print("from list 1: TrajectoryIterator")
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        xyz_list = traj.tolist()
        assert len(xyz_list) == traj.n_frames
        farray_0 = Trajectory(xyz_list, traj.top)
        for f0, f1 in izip(farray_0, traj):
            assert_almost_equal(f0.coords, f1.coords)

    @Timer()
    def test_8(self):
        #print("from list 2: Trajectory")
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        _farray = traj[:]
        xyz_list = _farray.tolist()
        assert len(xyz_list) == traj.n_frames
        farray_0 = Trajectory(xyz_list, traj.top)
        assert farray_0.size == traj.size
        for f0, f1 in izip(farray_0, traj):
            assert_almost_equal(f0.coords, f1.coords)

    @no_test
    @Timer()
    def test_9(self):
        #print("from list 3: Coord dataset")
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist = DataSetList()
        dslist.add_set("coords", "myname_crd")
        dslist.add_set("traj", "myname_traj")
        dslist[0].top = traj.top
        dslist[0].load(traj)

        dslist[1].top = traj.top.copy()
        dslist[1].load("./data/md1_prod.Tc5b.x")

        from pytraj.misc import get_atts
        #print(get_atts(dslist[0]))
        xyz_list = dslist[0].tolist()
        assert len(xyz_list) == traj.n_frames
        farray_0 = Trajectory(xyz_list, traj.top)
        for f0, f1 in izip(farray_0, traj):
            assert_almost_equal(f0.coords, f1.coords)

        xyz_list = dslist[1].to_ndarray()
        t = dslist[1]
        #print(t[0])
        xyz_list = dslist[1].xyz
        assert len(xyz_list) == traj.n_frames
        #print(farray_0)
        farray_0 = Trajectory()
        farray_0.top = traj.top.copy()
        farray_0.append_xyz(xyz_list)

        for f0, f1 in izip(farray_0, traj):
            assert_almost_equal(f0.coords, f1.coords)

    @test_if_having("mdtraj")
    def test_10(self):
        # from mdtraj
        import mdtraj as md
        from pytraj import Trajectory
        mtop = md.load_prmtop("./data/Tc5b.top")
        m_traj = md.load_mdcrd("./data/md1_prod.Tc5b.x", top=mtop)
        fa = Trajectory(m_traj.xyz, top="./data/Tc5b.top")
        aa_eq(fa.xyz, m_traj.xyz)

    def test_11(self):
        # from_iterable
        t = pt.Trajectory.from_iterable(TRAJ)
        aa_eq(t.xyz, TRAJ.xyz)

        t = pt.Trajectory.from_iterable(TRAJ(mask='@CA'))
        aa_eq(t.xyz, TRAJ['@CA'].xyz)

        t = pt.Trajectory.from_iterable(TRAJ(mask='@CA', rmsfit=(0, '@CA')))
        aa_eq(t.xyz, pt.get_coordinates(TRAJ(mask='@CA', rmsfit=(0, '@CA'))))

    def test_12_xyz(self):
        traj = TRAJ
        aa_eq(pt.Trajectory(xyz=traj.xyz,
                            top=traj.top.filename).xyz,
              pt.Trajectory(xyz=traj.xyz,
                            top=traj.top).xyz)


if __name__ == "__main__":
    unittest.main()
