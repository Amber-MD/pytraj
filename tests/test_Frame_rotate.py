from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca

# TODO : add assert


class Test(unittest.TestCase):
    def test_0(self):
        import numpy as np
        from pytraj.math import Matrix_3x3
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        f0 = traj[0].copy()
        frame = traj[0]
        mat = Matrix_3x3(range(9))
        npmat = mat.to_ndmatrix()

        # use both pytraj and numpy matrix
        frame.rotate_with_matrix(mat)
        f0.rotate_with_matrix(mat.as_ndmatrix())
        eq(f0.coords, frame.coords)


if __name__ == "__main__":
    unittest.main()
