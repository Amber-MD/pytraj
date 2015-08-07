from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        # TrajectoryIterator
        # assert will fail since we only allow to load a single
        # frame at a time to memory. Use frame_iter instead.
        aa_eq(pt.calc_molsurf([f for f in traj],
                              top=traj.top), pt.calc_molsurf(traj))

        # frame_iter
        aa_eq(pt.calc_molsurf([f for f in traj()],
                              top=traj.top), pt.calc_molsurf(traj))

        # chunk_iter: need to explicitly copy
        # aa_eq(pt.calc_molsurf([f for f in traj.chunk_iter(3)], top=traj.top),
        #       pt.calc_molsurf(traj))

        aa_eq(pt.calc_molsurf([f.copy() for f in traj.chunk_iter(3)],
                              top=traj.top), pt.calc_molsurf(traj))

        # split_iterators
        aa_eq(pt.calc_molsurf([f for f in traj.split_iterators(3)],
                              top=traj.top), pt.calc_molsurf(traj))

        # Trajectory
        aa_eq(pt.calc_molsurf([f for f in traj.to_mutable_trajectory()],
                              top=traj.top), pt.calc_molsurf(traj))


if __name__ == "__main__":
    unittest.main()
