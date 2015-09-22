from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.base import *
from pytraj import io
from pytraj.utils import eq, aa_eq
from pytraj.testing import cpptraj_test_dir, duplicate_traj
from pytraj.utils import goto_temp_folder


class TestWriteTraj(unittest.TestCase):
    def setUp(self):
        self.traj = pt.load_sample_data('tz2')

    def test_regular(self):
        traj = self.traj.copy()
        assert traj[0].has_box() == True

        with goto_temp_folder():
            # write traj with nobox info
            fname = "traj_nobox.nc"
            pt.write_traj(fname, traj, mode='nobox')
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
        fname = './output/test_xyz.nc'
        pt.write_traj(fname, xyz, top=self.traj.top, overwrite=True)
        t0 = pt.iterload(fname, top=self.traj.top)
        aa_eq(self.traj.xyz, t0.xyz)


if __name__ == "__main__":
    unittest.main()
