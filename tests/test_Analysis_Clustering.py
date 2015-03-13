import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.analyses.Analysis_Clustering import Analysis_Clustering

class Test(unittest.TestCase):
    def test_0(self):
        # TODO: Got error Error: cluster data set _DEFAULTCRD_ does not contain data.
        # FIXME: I know I did not read cpptraj code carefully :D
        top = mdio.load("./data/tz2.parm7")
        arg = """
        parm ./data/tz2.parm7
        trajin ./data/tz2.nc 
        #crdset C1
        cluster C1 :2-10 clusters 3 epsilon 4.0 summary ./output/avg.summary.dat nofit
        """
        act = Analysis_Clustering()
        act(arg, top)

if __name__ == "__main__":
    unittest.main()
