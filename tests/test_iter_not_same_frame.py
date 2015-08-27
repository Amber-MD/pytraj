from __future__ import print_function
import unittest
import pytraj as pt
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
        aa_eq(pt.molsurf([f for f in traj], top=traj.top), pt.molsurf(traj))

        # frame_iter
        aa_eq(pt.molsurf([f for f in traj()], top=traj.top), pt.molsurf(traj))

        # chunk_iter: need to explicitly copy
        # aa_eq(pt.molsurf([f for f in traj.iterchunk(3)], top=traj.top),
        #       pt.molsurf(traj))

        aa_eq(pt.molsurf([f.copy() for f in traj.iterchunk(3)],
                         top=traj.top), pt.molsurf(traj))

        # split_iterators
        aa_eq(pt.molsurf([f for f in traj.split_iterators(3)],
                         top=traj.top), pt.molsurf(traj))

        # Trajectory
        #print(pt.molsurf(traj))
        #print(pt.molsurf([f for f in traj], top=traj.top))
        #t0 = traj[:]
        #print(pt.molsurf([f.copy() for f in t0], top=traj.top))
        #print(pt.molsurf(traj[:]))

        # this does not work since the memory is freed to soon.
        # make a dummy object to traj?
        #aa_eq(pt.molsurf([f for f in traj[:]],
        #                  top=traj.top), pt.molsurf(traj))


if __name__ == "__main__":
    unittest.main()
