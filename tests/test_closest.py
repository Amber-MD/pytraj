#!/usr/bin/env python

from __future__ import print_function
import unittest
import numpy as np
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestClosest(unittest.TestCase):

    def test_closest(self):
        # raise if not has solvent
        traj0 = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")
        self.assertRaises(RuntimeError, lambda: pt.closest(traj0))

        traj = pt.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
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

        pt.write_traj('output/fi.nc', fi, top=top, overwrite=True)
        # load back
        t1 = pt.load('output/fi.nc', top=top)
        aa_eq(t0.xyz, t1.xyz)

        # make sure n_sovent=10 (default)
        n_solvents = 0
        for mol in top.mols:
            if mol.is_solvent():
                n_solvents += 1
        assert n_solvents == 10, 'must be 10 solvents'

        fi, top = pt.closest(traj)
        pt.write_traj('output/test.pdb', next(fi), top=top, overwrite=True)

    def test_closest_compared_to_cpptraj(self):
        traj = pt.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        state = pt.load_cpptraj_state('''
        autoimage
        closest 100 :1-13
        createcrd mycrd''', traj)
        state.run()

        fi, top = pt.closest(traj(autoimage=True), mask=':1-13', n_solvents=100)
        xyz = pt.get_coordinates(fi)
        t0 = pt.Trajectory(xyz=xyz, top=top)
        aa_eq(state.data['mycrd'].xyz, t0.xyz)

        # dtype = 'trajectory'
        t1 = pt.closest(traj(autoimage=True), mask=':1-13', n_solvents=100, dtype='trajectory')
        aa_eq(state.data['mycrd'].xyz, t1.xyz)


if __name__ == "__main__":
    unittest.main()
