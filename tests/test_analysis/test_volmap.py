from __future__ import absolute_import

import unittest
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq
from pytraj.all_actions import volmap

cm = "dummy.dat 0.5 0.5 0.5 :WAT@O buffer 2.0 centermask !:1-13 radscale 1.36"

txt = """
parm {}
trajin {} 1 1
rms first :1-13
center :1-13 mass origin
volmap {} {} {}
"""


class TestVolmap(unittest.TestCase):
    def test_volmap(self):
        traj = pt.iterload(fn('tz2.ortho.nc'), fn('tz2.ortho.parm7'))[:1]
        size = ''
        center = ''
        state = pt.load_cpptraj_state(
            txt.format(
                fn('tz2.ortho.parm7'), fn('tz2.ortho.nc'), cm, size, center))
        state.run()
        cpp_data = state.data[-1].values

        traj = traj.superpose(mask=':1-13').center(':1-13 mass origin')
        ds = pt.volmap(
            traj,
            mask=':WAT@O',
            grid_spacing=(0.5, 0.5, 0.5),
            buffer=2.0,
            centermask='!:1-13',
            radscale=1.36)

        aa_eq(cpp_data, ds)

        # assert
        self.assertRaises(
            AssertionError,
            lambda: pt.volmap(traj, mask=':WAT@O', grid_spacing='0.5 0.5 0.5'))
        self.assertRaises(
            AssertionError,
            lambda: pt.volmap(traj, mask=':WAT@O', grid_spacing=(0.5, 0.5)))
        self.assertRaises(ValueError, lambda: pt.volmap(traj, mask=':WAT@O', grid_spacing=(0.5, 0.5, 0.5), size='20 20 20'))

        # test size

    @unittest.skip('skip due to segmentation fault with cpptraj.OMP')
    def test_volmap_size(self):
        cm_no_buffer = cm.replace('buffer 2.0', '')
        state = pt.load_cpptraj_state(
            txt.format(
                fn('tz2.ortho.parm7'),
                fn('tz2.ortho.nc'), cm_no_buffer, 'size 20,20,20', ''))
        state.run()
        cpp_data = state.data[-1].values
        ds = volmap(
            traj,
            mask=':WAT@O',
            grid_spacing=(0.5, 0.5, 0.5),
            size=(20, 20, 20),
            buffer=2.0,
            centermask='!:1-13',
            radscale=1.36)
        aa_eq(cpp_data, ds)

    @unittest.skip('skip due to segmentation fault with cpptraj.OMP')
    def test_volmap_center(self):
        # test center
        cm_no_buffer = cm.replace('buffer 2.0', '')
        state = pt.load_cpptraj_state(
            txt.format(
                fn('tz2.ortho.parm7'),
                fn('tz2.ortho.nc'), cm_no_buffer, 'size 20,20,20',
                'center 0.5,0.5,0.5'))
        state.run()
        cpp_data = state.data[-1].values
        ds = volmap(
            traj,
            mask=':WAT@O',
            grid_spacing=(0.5, 0.5, 0.5),
            size=(20, 20, 20),
            center=(0.5, 0.5, 0.5),
            buffer=2.0,
            centermask='!:1-13',
            radscale=1.36)
        aa_eq(cpp_data, ds)

    def test_volmap_RuntimeError(self):
        # raise RuntimeError
        dry_traj = pt.iterload(fn('tz2.nc'), fn('tz2.parm7'))
        self.assertRaises(RuntimeError, lambda: pt.volmap(dry_traj, mask=':WAT@O',
                                                          grid_spacing=(0.5, 0.5, 0.5)))


if __name__ == "__main__":
    unittest.main()
