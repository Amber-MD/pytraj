from time import time
import unittest
from pytraj.base import *
from pytraj import TrajectoryIterator


class TestLoadingTime(unittest.TestCase):
    def test_0(self):
        def get_time(indices):
            t0 = time()
            traj = Trajectory(filename="./data/md1_prod.Tc5b.x",
                              top="./data/Tc5b.top",
                              indices=indices)
            return time() - t0

        def get_time_2(indices):
            traj = TrajectoryIterator(filename="./data/md1_prod.Tc5b.x",
                                      top="./data/Tc5b.top")
            t0 = time()
            traj[indices]
            return time() - t0

        def get_time_3(indices):
            traj = TrajectoryIterator(filename="./data/md1_prod.Tc5b.x",
                                      top="./data/Tc5b.top")
            t0 = time()
            traj[indices]
            return time() - t0

        #print("slice: ", get_time_2(slice(0, 5, 1)))
        #print(get_time(list(range(5))))
        #print(get_time(list(range(0, 5, 1))))
        #print("slice: ", get_time(slice(0, 5, 1)))
        #print(get_time(list(range(1)) + list(range(5, 2, -2))))
        #print(get_time(sorted(list(range(1)) + list(range(5, 2, -2)))))

        #print(get_time_3(slice(0, 5, 1)))

        traj = TrajectoryIterator(filename="./data/md1_prod.Tc5b.x",
                                  top="./data/Tc5b.top")
        #print(traj[:1])
        #print(traj[3:5:2])
        #print(traj[6:])
        ##print(traj[:1], traj[3:5:2], traj[6:])

        traj2 = Trajectory()
        traj2.top = traj.top.copy()
        traj2.join(traj[:1])
        traj2.join(traj[2:5])
        #print(traj2)


if __name__ == "__main__":
    unittest.main()
