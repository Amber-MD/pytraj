from __future__ import print_function
import unittest
from pytraj import DatasetList, analdict
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal as aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
from pytraj.datasets.DataSetList import DataSetList

class Test(unittest.TestCase):
    @test_if_having("numpy")
    def test_0(self):
        import numpy as np
        import pytraj.common_actions as pyca
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        # make sure not getting segmentation fault
        v0 = pyca.calc_vector(traj, "@CA @N,C,O")
        v1 = traj.calc_vector("@CA @N,C,O")

        # timecorr
        # http://ambermd.org/doc12/Amber15.pdf, page #619
        vec2 = pyca.calc_vector(traj, "v2 @18,@19,@20 corrplane", dtype='vector')
        act = analdict['timecorr']
        dslist = DataSetList()
        dslist.add_existing_set(vec2)

        # don't let python free memory for this `dslist`
        dslist.set_py_free_mem(False)
        command = "vec1 %s tstep 1.0 tcorr 100.0 out v2.out order 2" % vec2.name
        act(command, dslist=dslist)
        saved_data = np.loadtxt("./data/tc5b.myvec.out", skiprows=1).transpose()[1] 
        aa_eq(saved_data, dslist[-1].to_ndarray())


if __name__ == "__main__":
    unittest.main()
