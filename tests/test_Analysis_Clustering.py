import unittest
from pytraj.decorators import no_test
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.analyses.CpptrajAnalyses import Analysis_Clustering
from pytraj import DataSetList
from pytraj.datasets.DataSetList import DataSetList

class Test(unittest.TestCase):
    def test_0(self):
        dslist = DataSetList()
        dflist = DataFileList()

        dslist.add_set("traj", "my_clustering", "")
        from pytraj.trajs.Trajin_Single import Trajin_Single
        traj = Trajin_Single("./data/tz2.nc", "./data/tz2.parm7")
        dslist[0].top = traj.top
        dslist[0].add_trajin(traj)
        command = """
        crdset my_clustering :2-10 clusters 3 epsilon 4.0 summary ./output/avg.summary.dat nofit
        """
        act = Analysis_Clustering()
        act(command=command, top=traj.top, dslist=dslist, dflist=dflist)
        dflist.write_all_datafiles()

    def test_1(self):
        import numpy as np
        print ("use common_actions")
        from pytraj.common_actions import do_clustering

        traj = mdio.iterload("./data/tz2.nc", "./data/tz2.parm7")
        command = """
        :2-10 clusters 3 epsilon 4.0 summary ./output/avg.summary.do_clustering.dat nofit
        """
        dslist = do_clustering(traj, command, traj.top)
        self.assertIsInstance(dslist, DataSetList)
        print (dslist.to_dict())

        dslist = do_clustering(traj, command, traj.top, dtype='ndarray')
        self.assertIsInstance(dslist, np.ndarray)

if __name__ == "__main__":
    unittest.main()
