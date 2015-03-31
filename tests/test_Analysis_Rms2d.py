import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.analysis_dict import AnalysisDict
from pytraj.datasets.DataSet_Coords_TRJ import DataSet_Coords_TRJ
from pytraj.datasets.DataSet_Coords_CRD import DataSet_Coords_CRD

anadict = AnalysisDict()

class Test(unittest.TestCase):
    def test_0(self):
        print ("test DataSet_Coords_TRJ")
        dslist = DataSetList()
        dflist = DataFileList()

        trajin = "./data/md1_prod.Tc5b.x"
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist.add_set("traj", "test_traj", "my_default_name")
        dslist[0].top = traj.top
        dslist[0].add_trajin(traj)
        coords_traj = dslist[0]
        print (coords_traj.size)
        act = anadict['rms2d']
        print (dslist[0].name)
        act("@CA crdset test_traj rmsout ./output/_test_2drms.dat", traj.top,
            dslist=dslist,
            dflist=dflist)
        print (act)
        print (dslist.size)
        dflist.write_all_datafiles()

    def test_1(self):
        print ("test DataSet_Coords_CRD")
        dslist = DataSetList()
        dflist = DataFileList()

        trajin = "./data/md1_prod.Tc5b.x"
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        dslist.add_set("coords", "test_traj", "my_default_name")
        dslist[0].top = traj.top
        # add frame to crd_traj
        for frame in traj:
            dslist[0].append(frame)

        crd_traj = dslist[0]
        assert crd_traj.size == traj.size
        act = anadict['rms2d']
        print (dslist[0].name)
        act("@CA crdset test_traj rmsout ./output/_test_2drms_CRDtest.dat", traj.top,
            dslist=dslist,
            dflist=dflist)
        print (act)
        print (dslist.size)
        print (dslist[1].size)
        print (dslist[1].mkind)
        print (dslist[1].name)
        print (dslist[1].dtype)
        assert (dslist[1].dtype == 'matrix_flt')
        print (len(dslist[1][:]))
        #dflist.write_all_datafiles()

if __name__ == "__main__":
    unittest.main()
