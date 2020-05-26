#!/usr/bin/env python

from __future__ import print_function
import unittest
import numpy as np
import pytraj as pt
from utils import fn, outputname
from pytraj.testing import aa_eq
import pytest


class TestClosest(unittest.TestCase):
    def test_closest(self):
        # raise if not has solvent
        traj0 = pt.iterload(
            fn('tz2.nc'), fn('tz2.parm7'), frame_slice=[(0, 2)])
        with pytest.raises(RuntimeError):
            pt.closest(traj0)

        traj = pt.iterload(
            fn('tz2.ortho.nc'), fn('tz2.ortho.parm7'), frame_slice=[(0, 2)])
        fi, top = pt.closest(traj)

        coords = []
        for frame in fi:
            coords.append(frame.xyz.copy())
            assert isinstance(frame, pt.Frame), 'must be Frame'

        # make a Trajectory
        fi, top = pt.closest(traj)
        xyz = pt.get_coordinates(fi)
        t0 = pt.Trajectory(xyz=xyz, top=top)
        aa_eq(np.array(coords), t0.xyz)

        # test write to disk
        fi, top = pt.closest(traj)
        new_traj = pt.Trajectory(top=top)
        for frame in fi:
            new_traj.append(frame)

        pt.write_traj(outputname('fi.nc'), new_traj, overwrite=True)
        # load back
        t1 = pt.load(outputname('fi.nc'), top=top)
        aa_eq(t0.xyz, t1.xyz)

        # make sure n_sovent=10 (default)
        n_solvents = 0
        for mol in top.mols:
            if mol.is_solvent():
                n_solvents += 1
        assert n_solvents == 10, 'must be 10 solvents'

    def test_closest_compared_to_cpptraj(self):
        trajin = fn('tz2.ortho.nc')
        parm = fn('tz2.ortho.parm7')
        traj = pt.iterload(
            trajin, parm, frame_slice=[
                (0, 2),
            ])
        state = pt.load_cpptraj_state('''
        parm {}
        trajin {} 1 2
        autoimage
        closest 100 :1-13
        createcrd mycrd'''.format(parm, trajin))
        state.run()

        fi, top = pt.closest(
            traj(autoimage=True), mask=':1-13', n_solvents=100)
        xyz = pt.get_coordinates(fi)
        t0 = pt.Trajectory(xyz=xyz, top=top)
        aa_eq(state.data['mycrd'].xyz, t0.xyz)

        # dtype = 'trajectory'
        t1 = pt.closest(
            traj(autoimage=True),
            mask=':1-13',
            n_solvents=100,
            dtype='trajectory')
        aa_eq(state.data['mycrd'].xyz, t1.xyz)


if __name__ == "__main__":
    unittest.main()
