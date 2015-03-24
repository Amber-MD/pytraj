import unittest
from pytraj.base import *
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
        print(dset.data_format)
        print(dset.is_empty())
        print(dset.dtype)
        print(dset.column_width)
        db = cast_dataset(dset, dtype="general")
        #print cast_dataset.__doc__
        print(db.is_empty())
        #print dir(db)

    #@no_test
    def test_add_trajin(self):
        traj = DataSet_Coords_TRJ()
        traj.top = Topology("data/Tc5b.top")
        traj.load("data/md1_prod.Tc5b.x")

        trajin = TrajReadOnly(filename="data/md1_prod.Tc5b.x", top=traj.top)
        print(dir(trajin))
        dset = traj.alloc()
        db = cast_dataset(dset, dtype="general")
        print(db)

if __name__ == "__main__":
    unittest.main()

