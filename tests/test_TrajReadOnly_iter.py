import unittest
from pytraj.base import *


class TestTrajingIter(unittest.TestCase):
    def test_iter_0(self):
        from pytraj import TrajectoryIterator
        traj = TrajectoryIterator()
        traj.top = Topology("./data/Tc5b.top")
        traj.load(filename="data/md1_prod.Tc5b.x")

        assert traj.size == traj.n_frames == len(traj)

        for i in range(2):
            print()
            print("loop number = %s-th" % i)
            for frame in traj.frame_iter(start=0, stride=1):
                print(frame)
                assert frame.n_atoms == traj.top.n_atoms

        from time import time
        t0 = time()
        for frame in traj.frame_iter():
            print(frame)
        print(time() - t0)


if __name__ == '__main__':
    unittest.main()
