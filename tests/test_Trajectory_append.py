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

fa = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")[:]
traj = Trajectory(mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top"))

class Test(unittest.TestCase):
    def test_0(self):
        # test append
        traj = Trajectory()
        t = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        traj.top = t.top

        # append single Frame
        traj.append(t[0])
        assert traj.n_frames == 1

        # append xyz
        traj.append(t.xyz)
        assert traj.n_frames == t.n_frames + 1
        
        # append TrajReadOnly
        traj.append(t)
        assert traj.n_frames == t.n_frames * 2 + 1

        # append frame_iter
        traj.append(t.frame_iter())
        assert traj.n_frames == t.n_frames * 3 + 1

        # append _frame_iter_master
        from pytraj._shared_methods import _frame_iter_master
        traj.append(_frame_iter_master(t))
        assert traj.n_frames == t.n_frames * 4 + 1

        # append itself
        NFrames = traj.n_frames
        traj.append(traj)
        assert traj.n_frames == NFrames * 2

        # append itself frame_iter
        traj.append(traj.frame_iter(stop=2))
        assert traj.n_frames == NFrames * 2 + 3

        # append _frame_iter_master for itself
        NFrames = traj.n_frames
        traj.append(_frame_iter_master(traj))
        assert traj.n_frames == NFrames * 2

        # append _frame_iter_master for itself + other
        n0 = traj.n_frames
        n1 = t.n_frames
        traj.append(_frame_iter_master([traj, t]))
        assert traj.n_frames == 2 * n0 + n1
                

if __name__ == "__main__":
    unittest.main()
