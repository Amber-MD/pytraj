from __future__ import print_function
import unittest
from pytraj import *
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq


class Test(unittest.TestCase):
    def test_0(self):
        import numpy as np
        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))

        fa = traj[:]
        for i in range(3):
            fa.append_xyz(fa.xyz)

        xyz = fa.xyz / 10.
        fa.xyz = xyz
        aa_eq(xyz, fa.xyz)

        # try to build Trajectory from scratch
        fa2 = Trajectory()
        fa2.top = fa.top
        fa2._allocate(fa.n_frames, fa.n_atoms)
        fa2.xyz = fa.xyz[:]
        aa_eq(fa2.xyz, fa.xyz)

        # try to build Trajectory from scratch
        fa2 = Trajectory(xyz=fa.xyz, top=fa.top)
        aa_eq(fa2.xyz, fa.xyz)

        # timing
        np.empty((fa.n_frames, fa.n_atoms, 3))


if __name__ == "__main__":
    unittest.main()
