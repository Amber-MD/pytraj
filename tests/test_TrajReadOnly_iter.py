import unittest
from pytraj.base import *

max_frame = 999
class TestTrajingIter(unittest.TestCase):
    def test_iter_0(self):
        traj = TrajReadOnly()
        traj.load(filename="data/md1_prod.Tc5b.x", top=Topology("./data/Tc5b.top"))

        assert traj.size == traj.n_frames == len(traj)
        
        for i in range(2):
            print() 
            print("loop number = %s-th" % i)
            for frame in traj.frame_iter(start=0, chunk=1):
                assert frame.n_atoms == traj.top.n_atoms

        for frame in traj.frame_iter():
            print (frame)
            
if __name__ == '__main__':
    unittest.main()
