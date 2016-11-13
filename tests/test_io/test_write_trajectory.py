from __future__ import print_function
import unittest
from glob import glob
import pytraj as pt
from pytraj.testing import aa_eq
from pytraj.utils import tempfolder

from utils import fn

class TestWriteTraj(unittest.TestCase):

    def setUp(self):
        self.traj = pt.load_sample_data('tz2')

    def test_regular(self):
        traj = self.traj.copy()
        assert traj[0].has_box() == True

        with tempfolder():
            # write traj with nobox info
            fname = "traj_nobox.nc"
            pt.write_traj(fname, traj, options='nobox')
            t = pt.load(fname, traj.top)

            assert t[0].has_box() == False
            aa_eq(t.xyz, traj.xyz)

            # write from frame_iter, need to provide top
            fname = "traj_frame_iter.nc"
            # raise ValueError if there is no Topology
            pt.write_traj(fname, traj(), top=traj.top, overwrite=True)
            t = pt.iterload(fname, traj.top)
            aa_eq(t.xyz, traj.xyz)

    def test_write_xyz(self):
        xyz = self.traj.xyz
        fname = 'test_xyz.nc'
        with tempfolder():
            pt.write_traj(fname, xyz, top=self.traj.top, overwrite=True)
            t0 = pt.iterload(fname, top=self.traj.top)
            aa_eq(self.traj.xyz, t0.xyz)

    def test_split_and_write_traj(self):
        traj = pt.iterload([fn('Tc5b.x'), fn('Tc5b.x')], fn('Tc5b.top'))
        # duplcate
        assert traj.n_frames == 20
        top = traj.top

        with tempfolder():
            # test TrajectoryIterator object
            pt.tools.split_and_write_traj(traj,
                                          n_chunks=4,
                                          root_name='trajiterx',
                                          overwrite=True)
            flist = sorted(glob("trajiterx*"))
            traj4 = pt.iterload(flist, top)
            aa_eq(traj4.xyz, traj.xyz)

            # dcd ext
            pt.tools.split_and_write_traj(traj,
                                          4,
                                          root_name='ts',
                                          ext='dcd',
                                          overwrite=True)
            flist = sorted(glob("ts.*.dcd"))
            traj4 = pt.iterload(flist, top)
            aa_eq(traj4.xyz, traj.xyz)

    def test_raise(self):
        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))

        # list
        self.assertRaises(
            NotImplementedError,
            lambda: pt.write_traj('output/test.pdb', [traj[0], traj[1]], top=traj.top, frame_indices=range(3), overwrite=True))

        # single Frame
        self.assertRaises(
            ValueError,
            lambda: pt.write_traj('output/test.pdb', traj[0], top=traj.top, frame_indices=range(3), overwrite=True))

        # single Frame, no Topology
        self.assertRaises(
            ValueError,
            lambda: pt.write_traj('output/test.pdb', traj[0], top=None, overwrite=True))


if __name__ == "__main__":
    unittest.main()
