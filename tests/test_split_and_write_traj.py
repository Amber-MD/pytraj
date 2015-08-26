from __future__ import print_function
import unittest; import pytraj as pt
import pytraj as pt
from pytraj.utils import eq, aa_eq, eq_coords
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca
from pytraj.testing import goto_temp_folder


class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        # duplcate
        traj.load(traj.filename)
        assert traj.n_frames == 20
        top = traj.top

        # test Trajectory object
        traj2 = traj.to_mutable_trajectory()
        with goto_temp_folder():
            traj2.split_and_write_traj(n_chunks=4, root_name='trajx')
            # reload
            from glob import glob
            flist = sorted(glob("trajx*"))
            #print(flist)
            traj3 = pt.iterload(flist, top)
            aa_eq(traj2.xyz, traj3.xyz)
            #print(traj2.rmsd())
            #print(traj3.rmsd())

        # test TrajectoryIterator object
        traj.split_and_write_traj(n_chunks=4, root_name='./output/trajiterx')
        flist = sorted(glob("./output/trajiterx*"))
        #print(flist)
        traj4 = pt.iterload(flist, top)
        eq_coords(traj4, traj)

        # dcd ext
        traj.split_and_write_traj(4, root_name='./output/ts', ext='dcd')
        flist = sorted(glob("./output/ts.*.dcd"))
        #print(flist)
        traj4 = pt.iterload(flist, top)
        eq_coords(traj4, traj)


if __name__ == "__main__":
    unittest.main()
