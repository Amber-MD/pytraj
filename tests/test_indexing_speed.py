import unittest
from time import time
from pytraj.base import *
from pytraj.TrajReadOnly import TrajReadOnly

def get_time(indices):
    traj = TrajReadOnly(filename="./data/md1_prod.Tc5b.x", top="./data/Tc5b.top")
    t0 = time()
    print(traj.size)
    traj[indices]
    return time()- t0

traj = TrajReadOnly(filename="./data/md1_prod.Tc5b.x", top="./data/Tc5b.top")
N = traj.size
#traj[slice(0, N-1, 3)]

#print get_time(slice(None, N-1, 1))
#print get_time(slice(0, None, 1))
#print get_time(slice(-1, None, 1))

class TestSpeed(unittest.TestCase):
    def test_0(self):
        #f0 = traj[-3:-1][0]
        farray = traj[:-3:-1]
        f0 = farray[0]
        #print f0
        #f1 = traj[998]
        #print f1
        #
        traj2 = FrameArray(filename="./data/md1_prod.Tc5b.x", top="./data/Tc5b.top")
        f0 = traj2[-3:-1][0]
        print(f0)
        #f1 = traj2[998]
        #print f1
        
        farray2 = traj[:-3:-1] + traj[:] + traj[::-1]
        #farray2 += farray2
        print(farray2)

if __name__ == "__main__":
    unittest.main()
