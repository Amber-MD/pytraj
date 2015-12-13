#!/usr/bin/env python
from __future__ import print_function
import unittest
import os
import numpy as np
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.testing import cpptraj_test_dir


class TestGrid(unittest.TestCase):

    def test_0(self):
        from pytraj.math import Grid
        import numpy as np
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
        from pytraj.all_actions import calc_grid
        traj = pt.load_sample_data("tz2")[:]
        traj.autoimage()
        traj.rmsfit(mask=':1-13')
        d = calc_grid(traj, " 20 0.5 20 0.5 20 0.5 :WAT@O")

        d = calc_grid(traj(), " 20 0.5 20 0.5 20 0.5 :WAT@O", top=traj.top)

    def test_action_bounds(self):
        # creat mutable trajectory
        traj = pt.load('data/tz2.ortho.nc', 'data/tz2.ortho.parm7')
        pt.autoimage(traj)
        pt.superpose(traj, ref=0, mask=':1-13&!@H=', mass=True)
        grid_data = pt._grid(traj, mask=':1-13', grid_spacing=[0.5, 0., 0.])

        text = '''
        parm data/tz2.ortho.parm7
        trajin data/tz2.ortho.nc
        autoimage
        rms first :1-13&!@H= mass
        bounds :1-13 dx .5 name MyGrid
        '''

        state = pt.load_cpptraj_state(text)
        state.run()
        cpp_grid = state.data['MyGrid'].values
        aa_eq(cpp_grid, grid_data)

    def test_just_run_state(self):
        txt = '''
        parm data/tz2.truncoct.parm7
        trajin data/tz2.truncoct.nc
        reference data/tz2.truncoct.nc [REF]
        autoimage triclinic
        grid nonortho.dx boxref [REF] 50 50 50 :WAT@O pdb output/test.pdb
        '''

        state = pt.load_cpptraj_state(txt)
        state.run()


if __name__ == "__main__":
    unittest.main()
