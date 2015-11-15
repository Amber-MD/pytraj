from __future__ import absolute_import

import unittest
import numpy as np
import pytraj as pt
from pytraj.testing import aa_eq
from pytraj.common_actions import calc_volmap as volmap

cm = "output/dummy.dat 0.5 0.5 0.5 :WAT@O buffer 2.0 centermask !:1-13 radscale 1.36"

txt = """
parm data/tz2.ortho.parm7
trajin data/tz2.ortho.nc 1 1
rms first :1-13
center :1-13 mass origin
volmap {0}
""".format(cm)


class TestVolmap(unittest.TestCase):

    def test_volmap(self):
        traj = pt.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")[:1]
        state = pt.load_cpptraj_state(txt)
        state.run()
        cpp_data = state.data[-1].values

        traj = traj.superpose(mask=':1-13').center(':1-13 mass origin')
        ds = volmap(traj,
                    mask=':WAT@O',
                    grid_spacing='0.5 0.5 0.5',
                    buffer=2.0,
                    centermask='!:1-13',
                    radscale=1.36)

        aa_eq(cpp_data, ds)


if __name__ == "__main__":
    unittest.main()
