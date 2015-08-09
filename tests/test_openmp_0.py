import unittest
from pytraj.base import *
from pytraj import set_world_silent
from pytraj.utils import Timer
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.datasets.DataSet_Coords_CRD import DataSet_Coords_CRD
from pytraj.analyses.CpptrajAnalyses import Analysis_Rms2d


class Test(unittest.TestCase):
    def test_1(self):
        # just need to install libcpptraj with openmp
        # that's it

        # export OMP_NUM_THREADS=1
        # python ./test_openmp_0.py
        # export OMP_NUM_THREADS=8
        # python ./test_openmp_0.py

        dslist = DataSetList()
        dflist = DataFileList()

        trajin = "./data/md1_prod.Tc5b.x"
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        dslist.add_set("coords", "test_traj", "my_default_name")
        dslist[0].top = traj.top
        for i in range(45):
            dslist[0].load(traj)
        act = Analysis_Rms2d()

        with Timer() as t:
            act("crdset test_traj rmsout ./output/_test_2drms_CRDtest.openmp.dat",
                traj.top,
                dslist=dslist,
                dflist=dflist)
        print("time gap = %s" % t.time_gap())
        dflist.write_all_datafiles()

        # make sure to reproduce cpptraj to avoif false-impression :D
        import numpy as np
        matout = dslist[-1].get_full_matrix()

        tmp = np.loadtxt("./data/test_openmp.Tc5b.n_threads_1.dat",
                         skiprows=1,
                         usecols=range(1, dslist[0].size + 1))
        cpp_save = tmp.flatten()
        assert_almost_equal(cpp_save, matout)


if __name__ == "__main__":
    unittest.main()
