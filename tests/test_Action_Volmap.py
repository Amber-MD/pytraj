from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir

txt = """
parm ./data/tz2.ortho.parm7
trajin ./data/tz2.ortho.nc 1 1
rms first :1-13
center :1-13 mass origin 
volmap volmap.dx 0.5 0.5 0.5 :WAT@O buffer 2.0 centermask !:1-13 \\
radscale 1.36
volmap volmap2.dx 0.5 0.5 0.5 :WAT@O buffer 2.0 centermask !:1-13 \\
radscale 1.36 peakcut 0.10 peakfile peaks.xyz
"""


class Test(unittest.TestCase):
    @test_if_having("numpy")
    def test_0(self):
        import numpy as np
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        with open("./output/volmap.in", 'w') as f:
            f.write(txt)
        state = mdio.load_cpptraj_file("./output/volmap.in")
        #print(state.toplist[0])
        state.run()
        dslist = state.datasetlist
        for d0 in dslist._base_dataset_iter():
            #print(d0.dtype, d0.name)

        d0 = dslist[0]
        #print(d0)
        arr = dslist[0]
        #print(arr.shape)
        #print(arr)
        #print(dslist.get_legends())
        d1 = dslist[1]
        #print(d1.shape)
        mynp = d1.to_ndarray().flatten()
        mylist = d1.tolist()
        myview = d1.data
        #print(myview[0, 0, 0])
        #print(mylist[0][0])

        #print(np.where(mynp > 0))

    def test_1(self):
        from pytraj.common_actions import calc_volmap
        traj = mdio.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        ds = calc_volmap(
            traj,
            "0.5 0.5 0.5 :WAT@O buffer 2.0 centermask !:1-13 radscale 1.36 peakcut 0.10 peakfile peaks.xyz")
        #print(ds)
        ##print (ds.tolist())


if __name__ == "__main__":
    unittest.main()
