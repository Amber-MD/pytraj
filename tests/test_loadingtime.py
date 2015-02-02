from time import time
import unittest
from pytraj.base import *
from pytraj.TrajReadOnly import TrajReadOnly

class TestLoadingTime(unittest.TestCase):
    def test_0(self):
        def get_time(indices):
            t0 = time()
            traj = FrameArray(filename="./data/md1_prod.Tc5b.x", top="./data/Tc5b.top",
                                  indices=indices)
            return time() - t0
        
        def get_time_2(indices):
            traj = TrajReadOnly(filename="./data/md1_prod.Tc5b.x", top="./data/Tc5b.top")
            t0 = time()
            traj[indices]
            return time()- t0
        
        def get_time_3(indices):
            traj = TrajReadOnly(filename="./data/md1_prod.Tc5b.x", top="./data/Tc5b.top"
                   )
            t0 = time()
            traj[indices]
            return time()- t0
        
        print("slice: ", get_time_2(slice(0, 250, 1)))
        print(get_time(list(range(250))))
        print(get_time(list(range(0, 250, 1))))
        print("slice: ", get_time(slice(0, 250, 1)))
        print(get_time(list(range(100)) + list(range(500, 200, -2))))
        print(get_time(sorted(list(range(100)) + list(range(500, 200, -2)))))
        
        print(get_time_3(slice(0, 250, 1)))
        
        traj = TrajReadOnly(filename="./data/md1_prod.Tc5b.x", top="./data/Tc5b.top")
        print(traj[:100])
        print(traj[300:500:2])
        print(traj[600:])
        print(traj[:100] + traj[300:500:2] + traj[600:])
        
        traj2 = FrameArray()
        traj2.top = traj.top.copy()
        traj2 += traj[:100]
        traj2 += traj[200:500]
        print(traj2)
        
if __name__ == "__main__":
    unittest.main()
