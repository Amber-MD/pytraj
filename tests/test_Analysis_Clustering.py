import unittest
from pytraj.decorators import no_test
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.analyses.Analysis_Clustering import Analysis_Clustering

class Test(unittest.TestCase):
    def test_0(self):
        dslist = DataSetList()
        dflist = DataFileList()

        dslist.add_set("traj", "my_clustering", "")
        traj = mdio.load("./data/tz2.nc", "./data/tz2.parm7")
        dslist[0].top = traj.top
        dslist[0].add_trajin(traj)
        command = """
        crdset my_clustering :2-10 clusters 3 epsilon 4.0 summary ./output/avg.summary.dat nofit
        """
        act = Analysis_Clustering()
        act(command, traj.top, dslist=dslist, dflist=dflist)
        dflist.write_all_datafiles()

    def test_1(self):
        print ("use common_actions")
        from pytraj.common_actions import do_clustering

        traj = mdio.load("./data/tz2.nc", "./data/tz2.parm7")
        command = """
        :2-10 clusters 3 epsilon 4.0 summary ./output/avg.summary.do_clustering.dat nofit
        """
        do_clustering(command, traj, traj.top)

if __name__ == "__main__":
    unittest.main()
