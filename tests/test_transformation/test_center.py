from __future__ import print_function
import unittest
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq

from utils import fn


class TestCenter(unittest.TestCase):
    def test_center(self):
        traj = pt.iterload(fn('tz2.ortho.nc'), fn('tz2.ortho.parm7'))
        pt.load(fn("tz2.center_mass.nc"), traj.top)

        fa = traj[:]
        fa2 = traj[:].copy()
        pt.center(fa, mask=':1', mass=True)
        aa_eq(fa.xyz, fa.xyz, decimal=5)

        # raise if center not in 'origin', 'box'
        self.assertRaises(ValueError, lambda: pt.center(fa, center='oh'))

        # center to point

        pt.center(fa, ':1', center=[0, 1, 2], mass=True)
        aa_eq(pt.center_of_mass(fa, ':1')[0], [0, 1, 2])

        pt.center(fa, ':1', center=[0, 1, 2], mass=False)
        aa_eq(pt.center_of_geometry(fa, ':1')[0], [0, 1, 2])

        fa2.center(':1', center=[0, 1, 2], mass=False)
        aa_eq(pt.center_of_geometry(fa2, ':1')[0], [0, 1, 2])

        # on_disk
        traj_on_disk = pt.iterload(fn('tz2.ortho.nc'), fn('tz2.ortho.parm7'))
        pt.center(traj_on_disk, ':1', center=[0, 1, 2], mass=True)
        aa_eq(pt.center_of_mass(traj_on_disk, ':1')[0], [0, 1, 2])


if __name__ == "__main__":
    unittest.main()
