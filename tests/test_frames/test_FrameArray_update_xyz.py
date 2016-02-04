from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq, Timer
from pytraj.testing import cpptraj_test_dir


class Test(unittest.TestCase):

    def test_0(self):
        import numpy as np
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")

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
        xyz0 = np.empty((fa.n_frames, fa.n_atoms, 3))


if __name__ == "__main__":
    unittest.main()
