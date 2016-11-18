from __future__ import absolute_import
import unittest

import pytraj as pt
from pytraj.utils import Timer
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.datasets.c_datasetlist import DatasetList
from pytraj.datafiles import DataFileList
from pytraj.analysis.c_analysis.c_analysis import Analysis_Rms2d

from utils import fn


class Test(unittest.TestCase):

    def test_1(self):
        # just need to install libcpptraj with openmp
        # that's it

        # export OMP_NUM_THREADS=1
        # python ./test_openmp_0.py
        # export OMP_NUM_THREADS=8
        # python ./test_openmp_0.py

        dslist = DatasetList()
        dflist = DataFileList()

        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))

        dslist.add("coords", "test_traj")
        dslist[0].top = traj.top
        for i in range(45):
            dslist[0].load(traj.filename)
        act = Analysis_Rms2d()

        with Timer() as t:
            act("crdset test_traj rmsout ./output/_test_2drms_CRDtest.openmp.dat",
                dslist=dslist,
                dflist=dflist)

        # make sure to reproduce cpptraj to avoif false-impression :D
        import numpy as np
        matout = dslist[-1].get_full_matrix()

        tmp = np.loadtxt("./data/test_openmp.Tc5b.n_threads_1.dat",
                         skiprows=1,
                         usecols=range(1, dslist[0].size + 1))
        cpp_save = tmp.flatten()
        # use decimal = 3 to mathc cpptraj's format here
        assert_almost_equal(cpp_save, matout, decimal=3)


if __name__ == "__main__":
    unittest.main()
