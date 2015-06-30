from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca
from pytraj.compat import zip


class Test(unittest.TestCase):

    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        fa = traj[:]
        default_mass = [1.0 for _ in range(traj.top.n_atoms)]

        for frame in fa:
            aa_eq(default_mass, frame.mass)

        fa.set_frame_mass()
        fa2 = fa.copy()
        # make sure frame mass is copied too
        for f0, f1 in zip(fa, fa2):
            aa_eq(f0.mass, f1.mass)

        # make sure we get correct mass
        aa_eq(f0.mass, traj.top.mass)
        print(f0.mass)

if __name__ == "__main__":
    unittest.main()
