import unittest
from time import time
from pytraj.base import *
from pytraj.TrajectoryIterator import TrajectoryIterator


def get_time(indices):
    traj = TrajectoryIterator(
        filename="./data/md1_prod.Tc5b.x",
        top="./data/Tc5b.top")
    t0 = time()
    traj[indices]
    return time() - t0


traj = TrajectoryIterator(
    filename="./data/md1_prod.Tc5b.x",
    top="./data/Tc5b.top")
N = traj.n_frames


class TestSpeed(unittest.TestCase):
    def test_0(self):
        #f0 = traj[-3:-1][0]
        farray = traj[:-3:-1]
        f0 = farray[0]
        #f1 = traj[998]
        #
        traj2 = Trajectory(
            filename="./data/md1_prod.Tc5b.x",
            top="./data/Tc5b.top")
        f0 = traj2[-3:-1][0]
        #f1 = traj2[998]

        farray2 = Trajectory(top=traj2.top)
        farray2.join([traj[:-3:-1], traj[:], traj[::-1]])
        #farray2 += farray2


if __name__ == "__main__":
    unittest.main()
