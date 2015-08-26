from __future__ import print_function
import unittest; import pytraj as pt
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        from pytraj import adict

        #print(adict.keys())
        #print(adict['rmsd'])
        #print(adict['radgyr'])
        #print(adict['matrix'])

        act = adict['matrix']
        #print(act)

        act("",
            current_frame=(
                traj, traj(1, 5, 1), traj.chunk_iter(chunksize=2)),
            top=traj.top)
        #print(act.n_frames)
        assert act.n_frames == 24

        act("@CA", (traj, traj(1, 5, 1), traj.frame_iter(stride=2)), traj.top)
        #print(act.n_frames)
        assert act.n_frames == 43

        act("@CA", traj.chunk_iter(), traj.top)
        assert act.n_frames == 53


if __name__ == "__main__":
    unittest.main()
