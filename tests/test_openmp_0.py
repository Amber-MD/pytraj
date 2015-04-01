import unittest
from pytraj.base import *
from pytraj.utils import Timer
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.analysis_dict import AnalysisDict
from pytraj.datasets.DataSet_Coords_TRJ import DataSet_Coords_TRJ
from pytraj.datasets.DataSet_Coords_CRD import DataSet_Coords_CRD

anadict = AnalysisDict()

class Test(unittest.TestCase):

    def test_1(self):
        print ("test DataSet_Coords_CRD")
        dslist = DataSetList()
        dflist = DataFileList()

        trajin = "./data/md1_prod.Tc5b.x"
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        dslist.add_set("coords", "test_traj", "my_default_name")
        dslist[0].top = traj.top
        for i in range(130):
            dslist[0].load(traj)
        print (dslist[0].size)
        act = anadict['rms2d']
        print (dslist[0].name)

        with Timer() as t:
            act("crdset test_traj rmsout ./output/_test_2drms_CRDtest.openmp.dat", traj.top,
                dslist=dslist,
                dflist=dflist)
            dflist.write_all_datafiles()
        print ("time gap = %s" % t.time_gap())

        # make sure to reproduce cpptraj to avoif false-impression :D
        import numpy as np
        matout = dslist[-1].get_full_matrix()

        tmp = np.loadtxt("./data/test_openmp.Tc5b.n_threads_1.dat", skiprows=1, usecols=range(1, dslist[0].size+1))
        cpp_save = tmp.flatten()
        print (cpp_save[:20])
        print (matout[:20])
        print (cpp_save.__len__())
        print (matout.__len__())
        assert_almost_equal(cpp_save, matout)

if __name__ == "__main__":
    unittest.main()
