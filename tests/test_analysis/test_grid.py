#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import aa_eq

from utils import fn


class TestGrid(unittest.TestCase):

    def test_0(self):
        from pytraj.math import Grid
        nx = ny = nz = 3
        g = Grid(nx, ny, nz)
        assert g.size == nx**3
        assert g.nx == g.ny == g.nz == nx

        value = 1000.
        g[0, 0, 0] = value
        assert g[0, 0, 0] == value
        assert g._element(0, 0, 0) == value


class TestGridAction(unittest.TestCase):

    def test_action_grid(self):
        from pytraj.all_actions import grid
        traj = pt.load_sample_data("tz2")[:]
        traj.autoimage()
        traj.rmsfit(mask=':1-13')
        d = grid(traj, " 20 0.5 20 0.5 20 0.5 :WAT@O")

        d = grid(traj(), " 20 0.5 20 0.5 20 0.5 :WAT@O", top=traj.top)

    def test_action_bounds(self):
        # creat mutable trajectory
        traj = pt.load(fn('tz2.ortho.nc'), fn('tz2.ortho.parm7'))
        pt.autoimage(traj)
        pt.superpose(traj, ref=0, mask=':1-13&!@H=', mass=True)
        grid_data = pt._grid(traj, mask=':1-13', grid_spacing=[0.5, 0., 0.])

        text = '''
        parm {}
        trajin {}
        autoimage
        rms first :1-13&!@H= mass
        bounds :1-13 dx .5 name MyGrid
        '''.format(fn('tz2.ortho.parm7'),
                   fn('tz2.ortho.nc'))

        state = pt.load_cpptraj_state(text)
        state.run()
        cpp_grid = state.data['MyGrid'].values
        aa_eq(cpp_grid, grid_data)

    def test_just_run_state(self):
        txt = '''
        parm {}
        trajin {}
        reference {} [REF]
        autoimage triclinic
        grid nonortho.dx boxref [REF] 50 50 50 :WAT@O pdb output/test.pdb
        '''.format(fn('tz2.truncoct.parm7'),
                   fn('tz2.truncoct.nc'),
                   fn('tz2.truncoct.nc'))

        state = pt.load_cpptraj_state(txt)
        state.run()


if __name__ == "__main__":
    unittest.main()
