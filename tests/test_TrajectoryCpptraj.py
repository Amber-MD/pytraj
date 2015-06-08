from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir, duplicate_traj
import pytraj.common_actions as pyca
from pytraj.compat import zip

class Test(unittest.TestCase):
    def test_0(self):
        from pytraj.trajs.TrajectoryCpptraj import TrajectoryCpptraj
        from pytraj.trajs.Trajin_Single import Trajin_Single
        traj = Trajin_Single("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        tc = TrajectoryCpptraj(traj.filename, traj.top)
        assert traj.n_frames == tc.n_frames
        assert traj.top.n_atoms == tc.top.n_atoms

        assert TrajectoryCpptraj(top=traj.top.filename).n_atoms == traj.n_atoms
        assert TrajectoryCpptraj(top=traj.top).n_atoms == traj.n_atoms

        # __iter__
        for f0, f1 in zip(traj, tc):
            aa_eq(f0.xyz, f1.xyz)

        # __getitem__
        aa_eq(tc[0].xyz, traj[0].xyz)
        aa_eq(tc['@CA', 0].xyz, traj['@CA', 0].xyz)

        # frame_iter
        for f0, f1 in zip(traj(), tc()):
            aa_eq(f0.xyz, f1.xyz)

        # __len__
        assert len(tc) == tc.size == tc.n_frames

        # slice
        from pytraj.utils import Timer
        s = slice(0, 9, 2)
        @Timer()
        def normal_slice(tc):
            tc[s]

        @Timer()
        def fast_slice(tc):
            tc._fast_slice(s)

        normal_slice(tc)
        fast_slice(tc)

        tc.load_new("./data/md1_prod.Tc5b.x")
        aa_eq(tc.xyz, traj.xyz)
        tc.load("./data/md1_prod.Tc5b.x")
        assert tc.n_frames == 2 * traj.n_frames
        aa_eq(tc[:traj.n_frames].xyz, tc[traj.n_frames:].xyz)


if __name__ == "__main__":
    unittest.main()
