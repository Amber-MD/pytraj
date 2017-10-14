#!/usr/bin/env python

from __future__ import print_function
import unittest
from pytraj.testing import aa_eq
import pytraj as pt
from utils import fn


class TestImage(unittest.TestCase):
    def test_image(self):
        traj_on_disk = pt.iterload(fn('tz2.ortho.nc'), fn('tz2.ortho.parm7'))
        traj = traj_on_disk[:]

        command = 'origin center :WAT'

        traj = pt.image(traj, command)

        traj_on_disk.image(command)

        state = pt.load_cpptraj_state("""
        parm {}
        trajin {}
        image origin center :WAT
        createcrd mycrd
        """.format(fn('tz2.ortho.parm7'), fn('tz2.ortho.nc')))

        state.run()
        cpptraj_xyz = state.data['mycrd'].xyz

        aa_eq(cpptraj_xyz, traj.xyz)
        aa_eq(cpptraj_xyz, traj_on_disk.xyz)


if __name__ == "__main__":
    unittest.main()
