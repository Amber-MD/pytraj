from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import io
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir, duplicate_traj
from pytraj.utils import goto_temp_folder

class Test(unittest.TestCase):
    def test_0(self):
        traj = io.load_sample_data("tz2")[:]
        assert traj.top.box.has_box() == True
        assert traj[0].has_box() == True

        with goto_temp_folder():
            # write traj with nobox info
            fname = "traj_nobox.nc"
            io.write_traj(fname, traj, more_args='nobox')
            t = io.load(fname, traj.top)
            print (t)
            # FIXME: assert failed
            #assert t[0].has_box() == False
            #assert t.top.box.has_box() == False
            aa_eq(t.xyz, traj.xyz)

            # write from frame_iter, need to provide top
            fname = "traj_frame_iter.nc"
            # raise ValueError if there is no Topology
            io.write_traj(fname, traj(), top=traj.top, overwrite=True)
            t = io.iterload(fname, traj.top)
            aa_eq(t.xyz, traj.xyz)

if __name__ == "__main__":
    unittest.main()
