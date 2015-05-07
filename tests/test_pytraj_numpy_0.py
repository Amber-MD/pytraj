from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca
from pytraj.api import Trajectory
from pytraj.six_2 import izip

fa = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")[:]
traj = Trajectory(mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top"))

class Test(unittest.TestCase):
    def test_0(self):
        # test loading
        traj = Trajectory(mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top"))
        assert (traj.shape == traj.xyz.shape) 
        assert traj.ndim == traj.xyz.ndim

        for i, frame in enumerate(traj.frame_iter()):
            assert frame.n_atoms == traj.n_atoms == 304 # trp-cage
        assert i + 1 == traj.n_frames

    def test_1(self):
        # test append
        fa = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")[:]
        fa2 = fa.copy()
        fa2.join(fa)

        traj = Trajectory()
        traj.top = fa.top.copy()
        traj.append(fa)
        traj.append(fa)
        print (traj)
        assert traj.n_frames == fa.n_frames * 2
        assert traj.n_atoms == fa.n_atoms

        count = 0
        for f0, f1 in izip(fa2, traj):
            count += 1
            assert_almost_equal(f0.coords, f1.coords)
        assert count == traj.n_frames

    def test_calc(self):
        arr0 = fa.calc_radgyr().tolist()
        arr1 = traj.calc_radgyr().tolist()
        assert_almost_equal(arr0, arr1)

    def test_frame_iter(self):
        # frame_iter
        traj = Trajectory(mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")[:])
        for frame in traj.frame_iter(mask='@CA'):
            assert frame.n_atoms == 20

        traj2 = Trajectory()
        traj2.top = traj.top.copy()
        print (traj2)
        for frame in traj.frame_iter():
            print (frame)
            traj2.append(frame)
        #traj2.append(traj.frame_iter()) # infinite loop
        assert traj2.n_frames == traj.n_frames
        assert traj2.n_atoms == traj.n_atoms == 304

    @test_if_having("mdtraj")
    def test_3(self):
        # test mdtraj
        # TrajNumpy
        import mdtraj as md
        fa = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")[:]
        traj = Trajectory(fa)
        arr0 = md.rmsd(traj, traj, 0)
        f0 = fa[0].copy()

        arr1 = [frame.rmsd(f0) for frame in traj]
        assert_almost_equal(arr1, arr0)

    def test_4(self):
        # load
        traj = Trajectory()
        traj2 = Trajectory()
        fnames = ("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        trajread= mdio.iterload(*fnames)
        traj.top = trajread.top.copy()
        traj.load(fnames[0])
        traj2.append(traj)
        assert traj.n_frames == trajread.n_frames == traj2.n_frames == 10

        for f0, f1, f2 in izip(traj, traj2, trajread):
            assert f0.same_coords_as(f1) == True
            assert f0.same_coords_as(f2) == True

    def test_5(self):
        # test loading with filename + topology file
        traj = Trajectory("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        trajread = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        assert traj.n_frames == trajread.n_frames

        # test loading with filename + Topology object
        traj2 = Trajectory("./data/md1_prod.Tc5b.x", traj.top)
        assert traj2.n_frames == traj.n_frames

        # test loading with Topology object
        traj2 = Trajectory(top="./data/Tc5b.top")
        traj2.load("./data/md1_prod.Tc5b.x")
        assert traj2.n_frames == traj.n_frames

        # append
        old_n_frames = traj2.n_frames
        traj2.append(traj)
        assert traj2.n_frames == old_n_frames + traj.n_frames

        import numpy as np
        print (traj2[old_n_frames:].xyz.shape)
        #assert np.any(traj2[old_n_frames:].xyz.flatten(), traj.xyz.flatten()) == True

if __name__ == "__main__":
    unittest.main()
