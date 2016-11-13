#!/usr/bin/env python

from __future__ import print_function
import unittest
from pytraj.utils import aa_eq


class TestImage(unittest.TestCase):

    def test_image(self):
        import pytraj as pt
        traj_on_disk = pt.iterload('data/tz2.ortho.nc', 'data/tz2.ortho.parm7')
        traj = traj_on_disk[:]

        command = 'origin center :WAT'

        traj = pt.image(traj, command)

        traj_on_disk.image(command)

        state = pt.load_cpptraj_state("""
        parm data/tz2.ortho.parm7
        trajin data/tz2.ortho.nc
        image origin center :WAT
        createcrd mycrd
        """)

        state.run()
        cpptraj_xyz = state.data['mycrd'].xyz

        aa_eq(cpptraj_xyz, traj.xyz)
        aa_eq(cpptraj_xyz, traj_on_disk.xyz)

if __name__ == "__main__":
    unittest.main()
