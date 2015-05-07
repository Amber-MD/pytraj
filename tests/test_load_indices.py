import unittest
import numpy as np
from array import array
from pytraj.base import *
from pytraj.io import load
from pytraj.decorators import no_test
from pytraj.utils.check_and_assert import assert_almost_equal

from load_traj import load as npload

class TestIndices(unittest.TestCase):
    #@no_test
    def test_0(self):

        traj1 = TrajectoryIterator(filename="data/md1_prod.Tc5b.x", top="./data/Tc5b.top")
        print(traj1.size)
        indices = slice(9, 6, -1)
        print(indices.indices(traj1.size))
        traj0 = Trajectory()
        traj0.load(filename="./data/md1_prod.Tc5b.x", top=Topology("./data/Tc5b.top"), indices=indices)
        print(traj0.size)
        # load whole traj

        # check if loading correctly
        # if 'True': wrong indexing
        # we actually sorted indices before reading
        print(traj0[0][:10])
        # traj0[0] must has the same coords as [traj1[100]]
        assert traj0[0].same_coords_as(traj1[9]) == True
        assert traj0[1].same_coords_as(traj1[8]) == True
        assert traj0[2].same_coords_as(traj1[7]) == True

        print(traj0[0])
        assert traj0[0].rmsd(traj1[9]) < 1E-4

        rmsdlist = []
        ref = traj1[0].copy()
        for frame in traj1:
            #print frame
            rmsdlist.append(frame.rmsd(ref))

        nparr = np.array(rmsdlist)

        # make sure we don't suport other indices 
        traj2 = Trajectory()
        traj2.load(filename="./data/md1_prod.Tc5b.x", 
                   top=Topology("./data/Tc5b.top"), 
                   indices=list(range(4)) + list(range(9, 5, -1)) + [4,])
        assert traj2[-1].coords == traj1[4].coords


    def test_array_assigment(self):
        traj1 = TrajectoryIterator(filename="data/md1_prod.Tc5b.x", top="./data/Tc5b.top")[:]
        print(traj1[0][10])
        print(traj1[1][10])

        # assign traj1[0] 
        traj1[0] = traj1[1].copy()
        # make sure the assignment happed correctly
        assert traj1[0].same_coords_as(traj1[1]) == True

        print("update traj1[0] and make sure this does not affect traj[100]")
        traj1[0][10, 0] = 1000000.
        assert traj1[0][10, 0] == traj1[0, 10, 0] == 1000000.
        assert (traj1[0].same_coords_as(traj1[1])) == False
        assert traj1[0, 10, 0] != traj1[1, 10, 0]

        print(traj1[0][:11])
        print(traj1[1][:11])

    def test_1(self):
        traj0 = TrajectoryIterator(filename="data/md1_prod.Tc5b.x", top="./data/Tc5b.top")
        traj = TrajectoryIterator(filename="data/md1_prod.Tc5b.x", top="./data/Tc5b.top")[:]
        assert traj[0].coords == traj0[0].coords
        print(traj[0].coords[:10])

        traj2 = TrajectoryIterator(filename="data/md1_prod.Tc5b.x", top="./data/Tc5b.top")[:][:10]
        assert traj2[0].coords == traj0[0].coords

        traj.join(traj[:] + traj[0:100] + traj[9:3:-1])
        traj += traj[:]

        assert traj[0].coords != array('d', [0 for _ in range(traj[0].size)])
        assert traj[-1].coords != array('d', [0 for _ in range(traj[0].size)])

        for frame in traj:
            frame.zero_coords()

        assert traj[0].coords == array('d', [0 for _ in range(traj[0].size)])
        assert traj[-1].coords == array('d', [0 for _ in range(traj[0].size)])

    def test_del_top(self):
        # why here? lazy to make another file
        top = Topology("./data/Tc5b.top")
        top2 = top
        print(top2 == top)
        print(top2)
        print(top)
        del top

    def test_join_dummy(self):
        traj0 = TrajectoryIterator(filename="data/md1_prod.Tc5b.x", top="./data/Tc5b.top")[:]
        #traj0 += traj0
        traj0 += traj0[:]
        print(traj0)

    def test_load_indices_from_io(self):
        from pytraj import io as mdio
        traj0 = mdio.load(filename="data/md1_prod.Tc5b.x", top="./data/Tc5b.top", indices=(1, 3, 7))
        trajreadonly = mdio.iterload(filename="data/md1_prod.Tc5b.x", top="./data/Tc5b.top")

        assert isinstance(traj0, Trajectory)
        assert_almost_equal(traj0[0].coords, trajreadonly[1].coords)
        assert_almost_equal(traj0[1].coords, trajreadonly[3].coords)
        assert_almost_equal(traj0[2].coords, trajreadonly[7].coords)
        print (traj0)

if __name__ == "__main__":
    unittest.main()

