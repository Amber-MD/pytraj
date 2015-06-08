from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.datasets.DataSet_Coords_TRJ import DataSet_Coords_TRJ 
from pytraj.trajs.Trajin_Single import Trajin_Single

class Test(unittest.TestCase):
    def test_0(self):
        dslist = DataSetList()
        dslist0 = DataSetList()
        traj = Trajin_Single("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        act = adict['matrix']
        act("", traj, dslist=dslist)

        for ds in dslist:
            dslist0._add_copy_of_set(ds)

        assert dslist0.size == dslist.size

        dtraj = DataSet_Coords_TRJ()
        dtraj.top = traj.top.copy()
        dtraj.add_trajin(traj)
        print (dtraj)
        print ('dtraj.dtype = ', dtraj.dtype)
        assert dtraj.size == traj.size

        # add dtraj to dslist0
        dslist0._add_copy_of_set(dtraj)
        assert dslist0.size == 2

        # can we take the dtraj back?
        dtrajback = dslist0[1]

        # need to give topology
        dtrajback.top = dtraj.top
        print (dtrajback.size)

        for d in dslist0:
            print (d)

    def test_1(self):
        # FIXME: can not add coords_traj
        # I missed anything here?
        from pytraj import set_world_silent
        dslist = DataSetList()
        traj = Trajin_Single("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        coords_traj = DataSet_Coords_TRJ()
        coords_traj.top = traj.top
        coords_traj.add_trajin(traj)
        dslist.add_existing_set(coords_traj)
        dslist.add_existing_set(coords_traj)

    def test_2(self):
        dslist = DataSetList()
        dslist.add_set("traj", "my_name", "__my_default_name__")
        print (dslist.size)
        print (dslist[0].size)
        traj = Trajin_Single("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist[0].top = traj.top
        dslist[0].add_trajin(traj)
        assert dslist[0].size == traj.n_frames

if __name__ == "__main__":
    unittest.main()
