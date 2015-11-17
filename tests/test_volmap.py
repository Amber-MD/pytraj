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
volmap {0} {1}
"""


class TestVolmap(unittest.TestCase):

    def test_volmap(self):
        traj = pt.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")[:1]
        size = ''
        state = pt.load_cpptraj_state(txt.format(cm, size))
        state.run()
        cpp_data = state.data[-1].values

        traj = traj.superpose(mask=':1-13').center(':1-13 mass origin')
        ds = pt.volmap(traj,
                    mask=':WAT@O',
                    grid_spacing=(0.5, 0.5, 0.5),
                    buffer=2.0,
                    centermask='!:1-13',
                    radscale=1.36)

        aa_eq(cpp_data, ds)

        # assert
        self.assertRaises(AssertionError, lambda: pt.volmap(traj, mask=':WAT@O', grid_spacing='0.5 0.5 0.5'))
        self.assertRaises(AssertionError, lambda: pt.volmap(traj, mask=':WAT@O', grid_spacing=(0.5, 0.5)))
        self.assertRaises(ValueError, lambda: pt.volmap(traj, mask=':WAT@O', grid_spacing=(0.5, 0.5, 0.5), size='20 20 20'))

        # test size
        state = pt.load_cpptraj_state(txt.format(cm, 'size 20,20,20'))
        state.run()
        cpp_data = state.data[-1].values
        ds = volmap(traj,
                    mask=':WAT@O',
                    grid_spacing=(0.5, 0.5, 0.5),
                    size=(20, 20, 20),
                    buffer=2.0,
                    centermask='!:1-13',
                    radscale=1.36)
        aa_eq(cpp_data, ds)

        # raise RuntimeError
        dry_traj = pt.iterload('data/tz2.nc', 'data/tz2.parm7')
        self.assertRaises(RuntimeError, lambda: pt.volmap(dry_traj, mask=':WAT@O',
            grid_spacing=(0.5, 0.5, 0.5)))


if __name__ == "__main__":
    unittest.main()
