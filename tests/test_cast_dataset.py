import unittest
from pytraj.base import *
from pytraj.datasets import cast_dataset
from pytraj.decorators import no_test
from pytraj.datasets.DataSet_Coords_TRJ import DataSet_Coords_TRJ

class Test(unittest.TestCase):
    #@no_test
    def test_0(self):
        # DataSet_Coords_TRJ class
        traj = DataSet_Coords_TRJ()
        traj.top = Topology("data/Tc5b.top")
        traj.load("data/md1_prod.Tc5b.x")
        print(dir(traj))
        print(traj.size)
        dset = traj.alloc()
        #print dir(dset)
        print(dset.format)
        print(dset.is_empty())
        print(dset.dtype)
        print(dset.column_width)
        db = cast_dataset(dset, dtype="general")
        #print cast_dataset.__doc__
        print(db.is_empty())
        #print dir(db)

    #@no_test
    def test_add_trajin(self):
        dset_traj = DataSet_Coords_TRJ()
        dset_traj.top = Topology("data/Tc5b.top")
        dset_traj.load("data/md1_prod.Tc5b.x")

        # dummy casting, just want to make sure we get what we want
        db = cast_dataset(dset_traj, dtype="traj")
        db.top = dset_traj.top
        print(db)
        assert db.size == dset_traj.size
        print (db[0])

        # try to add to DataSetList
        dslist = DataSetList()
        dslist._add_copy_of_set(db)
        print (dslist[0])
        assert isinstance(dslist[0], DataSet_Coords_TRJ)

if __name__ == "__main__":
    unittest.main()

