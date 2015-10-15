import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
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
    def test_0(self):
        import numpy as np
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        with open("./output/volmap.in", 'w') as f:
            f.write(txt)
        state = mdio.load_cpptraj_file("./output/volmap.in")
        state.run()
        dslist = state.datasetlist
        for d0 in dslist._base_dataset_iter():
            pass

        d0 = dslist[0]
        arr = dslist[0]
        d1 = dslist[1]
        mynp = d1.to_ndarray().flatten()
        mylist = d1.tolist()
        myview = d1.data

    def test_1(self):
        from pytraj.common_actions import calc_volmap
        traj = mdio.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        ds = calc_volmap(
            traj,
            "0.5 0.5 0.5 :WAT@O buffer 2.0 centermask !:1-13 radscale 1.36 peakcut 0.10 peakfile peaks.xyz")


if __name__ == "__main__":
    unittest.main()
