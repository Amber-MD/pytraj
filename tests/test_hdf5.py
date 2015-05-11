import unittest
import numpy as np
from pytraj.base import *
from pytraj import io as io
from pytraj.utils import aa_eq
from pytraj.decorators import test_if_having

class Test(unittest.TestCase):
    @test_if_having("h5py")
    def test_0(self):
        traj = io.load_hdf5("./data/ala2.h5")
        print (traj)
        assert traj.top.has_box() == False

    @test_if_having("h5py")
    @test_if_having("mdtraj")
    def test_1(self):
        import mdtraj as md
        from mdtraj.testing import get_fn
        fn = get_fn("frame0.h5")
        m_traj = md.load(fn)
        print (m_traj)
        traj = io.load_hdf5(fn)
        traj2 = io.load_hdf5(fn, autoconvert=False)
        print (traj)
        assert traj.top.has_box() == True
        assert traj2.top.has_box() == True

        aa_eq(traj.xyz, m_traj.xyz * 10)
        aa_eq(traj2.xyz, m_traj.xyz)
        print (traj.top.box, traj2.top.box)
        save_list10 = [10.0, 10.0, 10.0, 90.0, 90.0, 90.0]
        save_list = [1.0, 1.0, 1.0, 90.0, 90.0, 90.0]
        blist = traj.top.box.tolist()
        assert blist == save_list10

        for frame in traj:
            blist = frame.box.tolist()
            assert blist == save_list10

        # no autoconvert
        for frame in traj2:
            blist = frame.box.tolist()
            assert blist == save_list

    @test_if_having("h5py")
    @test_if_having("mdtraj")
    def test_1(self):
        import mdtraj as md
        import numpy as np
        from pytraj import Trajectory
        from mdtraj.testing import get_fn
        from timeit import timeit
        fn = get_fn("frame0.h5")
        m_traj = md.load(fn)

        def mdtraj_load():
            md.load(fn)

        def pytraj_load():
            io.load_hdf5(fn)

        def pytraj_convert():
            io.load_mdtraj(m_traj)

        fa = Trajectory()
        n_frames, n_atoms, _ = m_traj.xyz.shape
        fa._allocate(n_frames, n_atoms)
        crd = m_traj.xyz.astype(np.float64)

        def _allocate():
            fa.update_xyz(crd)

        def assign_xyz():
            _xyz = np.empty_like(m_traj.xyz)
            _xyz[:] = crd

        def use_api():
            io.load_hdf5(fn, restype='api.Trajectory')

        t_mdtraj = timeit(mdtraj_load, number=10)
        t_pytraj = timeit(pytraj_load, number=10)
        t_convert = timeit(pytraj_convert, number=10)
        t_alloc = timeit(_allocate, number=10)
        t_xyz = timeit(assign_xyz, number=10)
        t_api = timeit(use_api, number=10)

        print (t_mdtraj, t_pytraj, t_convert, t_alloc, t_xyz, t_api)

    @test_if_having("h5py")
    def test_2(self):
        # test read from buffer
        import h5py
        traj = io.load_hdf5("./data/ala2.h5")
        with h5py.File("./data/ala2.h5") as fh:
            traj2 = io.load_hdf5(fh)
            aa_eq(traj2.xyz, traj.xyz)

if __name__ == "__main__":
    unittest.main()
