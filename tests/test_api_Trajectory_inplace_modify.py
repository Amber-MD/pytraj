from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        farray = traj[:]

        t0api = pt.api.Trajectory(traj)
        aa_eq(farray.unitcells, t0api.unitcells)

        # autoimage
        farray.autoimage()
        t0api.autoimage()
        aa_eq(farray.xyz, t0api.xyz)

        # rotate_dihedral
        pt.rotate_dihedral(t0api, '3:phi:120')
        pt.rotate_dihedral(farray, '3:phi:120')
        aa_eq(farray.xyz, t0api.xyz)

        aa_eq(pt.calc_phi(t0api, '3').values, [120
                                               for _ in range(t0api.n_frames)])


if __name__ == "__main__":
    unittest.main()
