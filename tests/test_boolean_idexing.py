from __future__ import print_function
import unittest; import pytraj as pt
import pytraj as pt
from pytraj.utils import eq, aa_eq, eq_coords
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        import numpy as np
        for _ in range(5):
            # get different random number
            # TrajectoryIterator (immutable): use `iterload`
            traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

            # make sure to raise if mixing dtype
            self.assertRaises(NotImplementedError, lambda: traj[[True, 2]])
            self.assertRaises(NotImplementedError, lambda: traj[:][[True, 2]])

            assert len(traj[[True, False]]) == 1
            assert pt.tools.rmsd(traj[[True, False]].xyz, traj[0].xyz,
                                 flatten=True) < 1E-6

            brr = np.random.randint(0, 2, traj.n_frames) > 0
            arr = np.arange(traj.n_frames)[brr]
            print(arr)
            assert pt.tools.rmsd(
                traj[arr].xyz, traj[brr].xyz,
                flatten=True) < 1E-6

            # Trajectory (mutable)
            traj = pt.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

            assert len(traj[[True, False]]) == 1
            assert pt.tools.rmsd(traj[[True, False]].xyz, traj[0].xyz,
                                 flatten=True) < 1E-6

            brr = np.random.randint(0, 2, traj.n_frames) > 0
            arr = np.arange(traj.n_frames)[brr]
            print(arr)
            assert pt.tools.rmsd(
                traj[arr].xyz, traj[brr].xyz,
                flatten=True) < 1E-6


if __name__ == "__main__":
    unittest.main()
