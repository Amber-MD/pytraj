import unittest
from pytraj.six_2 import izip
from pytraj.base import *
from pytraj.datasets.DataSet_Coords_CRD import DataSet_Coords_CRD as DataSet_Coords_CRD
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj import adict
import numpy as np

class Test(unittest.TestCase):
    def test_0(self):
        TRAJ0 = TrajReadOnly(filename="./data/md1_prod.Tc5b.x", top="./data/Tc5b.top")
        print(TRAJ0.size)
        traj = DataSet_Coords_CRD()
        traj.top = TRAJ0.top.copy()

        # adding frames to DataSet_Coords_CRD
        Nframe = 10
        for i in range(Nframe):
            traj.append(TRAJ0[i])
        assert traj.size == Nframe == traj.n_frames
        print (dir(traj))

        # test negative indexing
        assert_almost_equal(traj[-1].coords, traj[Nframe-1].coords)
        assert_almost_equal(traj[-2].coords, traj[Nframe-2].coords)

        # make sure to add correct frames from TRAJ0
        for f0, f1 in izip(TRAJ0, traj):
            print ("coords f0 f1")
            print (f0[0], f1[0])
            assert_almost_equal(f0.coords, f1.coords)

        # test iteration
        for f1 in traj:
            print (f1)

        # test set Frame
        traj[0] = traj[9]
        assert_almost_equal(traj[0].coords, traj[9].coords)
        print (traj[0])
        
    def test_1(self):
        "test Action"
        TRAJ0 = TrajReadOnly(filename="./data/md1_prod.Tc5b.x", top="./data/Tc5b.top")
        traj = DataSet_Coords_CRD()
        traj.top = TRAJ0.top.copy()

        Nframe = 10
        for i in range(Nframe):
            traj.append(TRAJ0[i])

        # Action testing
        dslist = DataSetList()

        act = adict['distance']
        act(":2@CA :10@CA", traj, TRAJ0.top, dslist)
        d0 = dslist[0]
        print (d0.size)
        print (d0[:])

        d1 = np.loadtxt("./data/CAres2_CAres10.Tc5b.dat", skiprows=1).transpose()[1]
        print (d1)
        assert_almost_equal(dslist[0][:], d1)

    def test_2(self):
        from pytraj import cast_dataset
        from pytraj.datasets.DataSet import DataSet
        from pytraj.datasets.DataSet_Coords_CRD import DataSet_Coords_CRD
        "test cast_dataset"
        TRAJ0 = TrajReadOnly(filename="./data/md1_prod.Tc5b.x", top="./data/Tc5b.top")
        traj = DataSet_Coords_CRD()
        traj.top = TRAJ0.top.copy()

        Nframe = 10
        for i in range(Nframe):
            traj.append(TRAJ0[i])

        dataset = traj.alloc()
        assert isinstance(dataset, DataSet) == True

        # cast back
        casted_traj = cast_dataset(dataset, 'coords')
        print (casted_traj)
        assert isinstance(casted_traj, DataSet_Coords_CRD) == True

        # FIXME: make sure we get frames correctly after casting
        print ("test more")
        print (casted_traj.size) # = 0, wrong
        print (dir(casted_traj))
        #assert casted_traj.size == TRAJ0.size

        #for f0, f1 in izip(TRAJ0, casted_traj):
        #    print (f0, f1)
        #    assert_almost_equal(f0.coords, f1.coords)
        
    def test_3(self):
        # test commond DataSet methods
        TRAJ0 = TrajReadOnly(filename="./data/md1_prod.Tc5b.x", top="./data/Tc5b.top")
        print(TRAJ0.size)
        traj = DataSet_Coords_CRD()
        traj.top = TRAJ0.top.copy()

        # adding frames to DataSet_Coords_CRD
        Nframe = 10
        for i in range(Nframe):
            traj.append(TRAJ0[i])
            f0 = TRAJ0[i]
            f1 = traj[-1]
            assert_almost_equal (f0.coords, f1.coords)
        assert traj.size == Nframe == traj.n_frames

        dset = traj
        print (dset.is_empty())
        fcp = traj[0].copy()
        assert f0.n_atoms == TRAJ0[0].n_atoms

        traj.append(fcp)
        assert traj.n_frames == TRAJ0.n_frames + 1
        print (traj[-1].coords[:10])
        print (traj[0].coords[:10])
        print (traj[-1][0])
        print (traj[0][0])
        assert_almost_equal (traj[-1].coords, traj[0].coords)

        assert traj.dtype == 'coords'

if __name__ == "__main__":
    unittest.main()
