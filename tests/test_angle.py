from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca
from pytraj.misc import from_legends_to_indices


class Test(unittest.TestCase):
    def test_0(self):
        import numpy as np
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        fa = traj[:]
        mask = ':1@CA :14@CB :15CA'
        d0 = pyca.calc_angle(traj, mask, dtype='dataset').to_ndarray()
        d1 = traj.calc_angle(mask)
        d2 = fa.calc_angle(mask)

        aa_eq(d0, d1)
        aa_eq(d0, d2)

        Nsize = 10
        arr = np.random.randint(0, 300, size=Nsize * 3).reshape(Nsize, 3)
        d3 = fa.calc_angle(arr)
        d4 = traj.calc_angle(arr)
        d5 = pyca.calc_angle(traj, arr)
        d6 = pyca.calc_angle(fa, arr)
        d7 = pyca.calc_angle([fa, traj], arr, n_frames=2 * fa.n_frames)
        aa_eq(d3, d4)
        aa_eq(d3, d5)
        aa_eq(d3, d6)
        aa_eq(d3.T, d7.T[:fa.n_frames])
        aa_eq(d3.T, d7.T[fa.n_frames:])
        print(d3)


if __name__ == "__main__":
    unittest.main()
