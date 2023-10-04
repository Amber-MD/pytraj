from __future__ import print_function
import unittest
from pytraj import *
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq
import numpy as np


class Test(unittest.TestCase):
    # tests to see if function can find an atom if we search with a coordinate
    # very close to the atom's xyz coordinates
    def test_0(self):
        tz2_traj = pt.datafiles.load_tz2()

        # find coordinates of a given atom, make sure count_in_voxel finds
        # that atom at the right frame
        idx = 0
        frame = 0
        xyz = tz2_traj[frame].atom(idx)
        pop = pt.count_in_voxel(tz2_traj, tz2_traj.top, "", xyz, 3)
        assert (idx in pop[frame])

        idx = 3
        frame = 4
        xyz = tz2_traj[frame].atom(idx)
        pop = pt.count_in_voxel(tz2_traj, tz2_traj.top, "", xyz, 3)
        assert (idx in pop[frame])

    # tests to make sure function doesn't put an atom in the list if its
    # not in the voxel
    def test_1(self):
        tz2_traj = pt.datafiles.load_tz2()
        idx = 0
        x, y, z = tz2_traj[0].atom(idx)
        size = 3
        voxel = (x + size, y + size, z + size)
        pop = pt.count_in_voxel(tz2_traj, tz2_traj.top, "", voxel, size)

        assert (idx not in pop[0])


if __name__ == "__main__":
    unittest.main()
