#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestImage(unittest.TestCase):

    def test_image(self):
        import pytraj as pt
        traj_on_disk = pt.iterload('data/tz2.ortho.nc', 'data/tz2.ortho.parm7')
        traj = traj_on_disk[:]

        traj = pt.image(traj, 'origin center :WAT')

        state = pt.load_cpptraj_state("""
        parm data/tz2.ortho.parm7
        trajin data/tz2.ortho.nc
        image origin center :WAT
        createcrd mycrd
        """)

        state.run()
        cpptraj_xyz = state.data['mycrd'].xyz

        aa_eq(cpptraj_xyz, traj.xyz)

if __name__ == "__main__":
    unittest.main()
