import unittest
from pytraj.base import *
from pytraj.datasets.DataSet_Coords_CRD import DataSet_Coords_CRD as Trajectory
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        TRAJ = TrajReadOnly(filename="./data/md1_prod.Tc5b.x", top="./data/Tc5b.top")
        print(TRAJ.size)
        traj = Trajectory()
        traj.top = TRAJ.top.copy()

        Nframe = 10
        for i in range(Nframe):
            traj.append(TRAJ[i])
        assert traj.size == Nframe

        # test negative indexing
        assert_almost_equal(traj[-1][0], traj[Nframe-1][0])
        assert_almost_equal(traj[-2][0], traj[Nframe-2][0])

        # test set Frame
        traj[0] = traj[9]
        assert_almost_equal(traj[0][0], traj[9][0])
        
if __name__ == "__main__":
    unittest.main()
