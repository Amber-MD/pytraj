from __future__ import print_function
from glob import glob
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq, eq_coords
import pytraj.common_actions as pyca
from pytraj.testing import goto_temp_folder


class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        # duplcate
        traj.load(traj.filename)
        assert traj.n_frames == 20
        top = traj.top

        # test TrajectoryIterator object
        pt.tools.split_and_write_traj(traj,
                                      n_chunks=4,
                                      root_name='./output/trajiterx')
        flist = sorted(glob("./output/trajiterx*"))
        traj4 = pt.iterload(flist, top)
        eq_coords(traj4, traj)

        # dcd ext
        pt.tools.split_and_write_traj(traj, 4,
                                      root_name='./output/ts',
                                      ext='dcd')
        flist = sorted(glob("./output/ts.*.dcd"))
        traj4 = pt.iterload(flist, top)
        eq_coords(traj4, traj)


if __name__ == "__main__":
    unittest.main()
