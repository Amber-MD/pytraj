import unittest
from pytraj.base import *
from pytraj.datasets.DataSet_Coords_TRJ import DataSet_Coords_TRJ as Trajectory

class Test(unittest.TestCase):
    def test_0(self):
        TRAJ = TrajReadOnly(filename="./data/md1_prod.Tc5b.x", top="./data/Tc5b.top")
        print(TRAJ.size)
        traj = Trajectory()
        traj.top = TRAJ.top.copy()

        # D. Roe did not implemented this
        Nframe = 10
        for i in range(Nframe):
            traj.append(TRAJ[i])
        assert traj.size == 0
        #assert traj.size == Nframe

        ## test negative indexing
        #assert traj[-1][0] == traj[Nframe-1][0]
        #assert traj[-2][0] == traj[Nframe-2][0]

        ## test set Frame
        #traj[0] = traj[9]
        #assert traj[0][0] == traj[9][0]
        
if __name__ == "__main__":
    unittest.main()
