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
            for farray in traj.frame_iter(start=0, chunk=9):
                pass
            
            print(farray)
            print(farray[0])
            print(farray[0].coords[:10])
            print(farray[0].buffer1d)
            print(farray[0][0])

if __name__ == '__main__':
    unittest.main()
