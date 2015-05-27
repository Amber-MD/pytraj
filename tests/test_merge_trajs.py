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
        from pytraj.testing import make_fake_traj
        from pytraj.misc import merge_trajs
        import numpy as np
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        traj0 = traj[:5].copy()
        traj1 = traj[5:].copy()
        print (traj0, traj1)
        traj_merged  =  merge_trajs(traj0, traj1, start_new_mol=True)
        print (traj0, traj1)
        print (traj_merged)
        assert traj_merged.shape  == (traj0.n_frames, traj0.n_atoms + traj1.n_atoms, 3)

        n_atoms0 = traj0.n_atoms
        n_atoms_final = traj_merged.n_atoms

        for f0, f1, frame in zip(traj0, traj1, traj_merged):
            aa_eq(f0.xyz, frame.xyz[:n_atoms0])
            aa_eq(f1.xyz, frame.xyz[n_atoms0:])

        assert traj_merged.top.n_atoms == 2 * n_atoms0
        #assert traj_merged.top.n_mols == 2 #

        traj_merged_2 = merge_trajs(traj0, traj1, start_new_mol=False)
        assert traj_merged_2.top.n_mols == 1

        # merge from frame_iter
        traj_merged = merge_trajs((traj(1, 5), traj.top),
                                  (traj(1, 10, 2), traj.top),
                                  n_frames = 5)
        print (traj0.top, traj1.top)
        saved_fiter_0 = mdio._load_from_frame_iter(traj(1, 5), top=traj.top)
        saved_fiter_1 = mdio._load_from_frame_iter(traj(1, 10, 2), top=traj.top)

        for f0, f1, frame in zip(saved_fiter_0, saved_fiter_1, traj_merged):
            aa_eq(f0.xyz, frame.xyz[:n_atoms0])
            aa_eq(f1.xyz, frame.xyz[n_atoms0:])


if __name__ == "__main__":
    unittest.main()
