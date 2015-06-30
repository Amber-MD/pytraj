from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca
from pytraj.utils import Timer
from pytraj.compat import range


class Test(unittest.TestCase):

    @test_if_having("h5py")
    @test_if_having("mdtraj")
    def test_0(self):
        from mdtraj.testing import get_fn
        import mdtraj
        from pytraj.io import load_hdf5
        import h5py
        traj = io.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        fa = traj[:]
        fa[""] = traj.xyz
        fn = get_fn("frame0.h5")
        m0 = mdtraj.load(fn)
        top = io.load_mdtraj(m0).top

        @Timer()
        def normal_load():
            io.load_hdf5(fn)

        @Timer()
        def load_from_Trajectory():
            Trajectory(m0, top)

        @Timer()
        def _tfancy_load():
            _tfancy = Trajectory(top=top)
            _tfancy._allocate(m0.n_frames, m0.n_atoms)
            _tfancy['*'] = m0.xyz.astype('f8')

        @Timer()
        def _mdtraj_load():
            m = mdtraj.load(fn)

        print("normal load")
        normal_load()
        print("load_from_Trajectory")
        load_from_Trajectory()
        print("_tfancy")
        _tfancy_load()
        print("_mdtraj_load")
        _mdtraj_load()

        _t = Trajectory(m0, top)
        _m = mdtraj.load(fn)
        aa_eq(_t.xyz, _m.xyz * 10)

        _tfancy = Trajectory(top=top)
        _tfancy._allocate(m0.n_frames, m0.n_atoms)
        _tfancy['*'] = m0.xyz.astype('f8')
        aa_eq(_tfancy.xyz, _m.xyz)

    @test_if_having("mdtraj")
    def test_0(self):
        # test load hdf5 with given Topology
        from mdtraj.testing import get_fn
        import mdtraj
        from pytraj.io import load_hdf5
        import h5py
        traj = io.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        fa = traj[:]
        fn = get_fn("frame0.h5")
        m0 = mdtraj.load(fn)
        ptop = io.load_mdtraj(m0).top

        @Timer()
        def test_load_hdf5_no_top(n_times):
            for i in range(n_times):
                io.load_hdf5(fn)

        @Timer()
        def test_load_hdf5_with_top(n_times, ptop):
            for i in range(n_times):
                io.load_hdf5(fn, top=ptop)

        n_times = 50
        print("test_load_hdf5_no_top")
        test_load_hdf5_no_top(n_times)
        print("test_load_hdf5_with_top")
        test_load_hdf5_with_top(n_times, ptop)
        print(ptop)

if __name__ == "__main__":
    unittest.main()
