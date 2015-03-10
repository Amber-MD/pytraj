import unittest
from pytraj.base import *
from pytraj.decorators import no_test
from pytraj.datasets.DataSet_Coords_TRJ import DataSet_Coords_TRJ as Trajectory
from pytraj.datasets.DataSet_Coords_TRJ import DataSet_Coords_TRJ
from pytraj.utils import assert_almost_equal

class Test(unittest.TestCase):
    @no_test
    def test_0(self):
        TRAJ = TrajReadOnly(filename="./data/md1_prod.Tc5b.x", top="./data/Tc5b.top")
        print(TRAJ.size)

        # creat DataSet_Coords_TRJ object
        # man, should we use FrameArray?
        # might not be a big matter since adding frames are not that expensive vs clustering
        traj = Trajectory()
        traj.top = TRAJ.top.copy()

        # DataSet_Coords_TRJ does not have GetFrame method
        # what should I do?
        Nframe = 10
        for frame in TRAJ:
            traj.append(frame)

        # seriously I need to set size for DataSet_Coords_TRJ?
        traj.size = TRAJ.size
        assert traj.size == Nframe
        for f0 in traj:
            print (f0)
        #print (traj[0])
        #frame = Frame()
        #print traj.get_frame(3, frame, traj.top)
        #print (frame)

        ## test negative indexing
        #assert traj[-1][0] == traj[Nframe-1][0]
        #assert traj[-2][0] == traj[Nframe-2][0]

        ## test set Frame
        #traj[0] = traj[9]
        #assert traj[0][0] == traj[9][0]

    def test_1(self):
        traj = DataSet_Coords_TRJ()
        traj.top = Topology("./data/Tc5b.top")
        traj.load("./data/md1_prod.Tc5b.x")
        assert traj.size == 10
        traj.load("./data/md1_prod.Tc5b.x", arg="1 5 2")
        assert traj.size == 13

        for f0 in traj:
            print (f0)

        print ("getitem")
        print (traj[12])

        print (traj)
        ds0 = traj.alloc()
        print (ds0)
        print (dir(ds0))

        # try perform action
        from pytraj.common_actions import calc_distance
        d0 = calc_distance(":2@CA :10@CA", traj)
        print (d0.size)
        print (d0[:])

        # make sure we load correct frames
        assert_almost_equal(traj[0].coords, traj[10].coords)
        assert d0[0] == d0[10]
        
if __name__ == "__main__":
    unittest.main()
