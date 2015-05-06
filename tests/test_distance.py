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
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        fa = traj[:]
        mask = ':1@CA :14@CB'
        d0 = pyca.calc_distance(traj, mask).to_ndarray()
        d1 = traj.calc_distance(mask).to_ndarray()
        d2 = fa.calc_distance(mask).to_ndarray()

        aa_eq(d0, d1)
        aa_eq(d0, d2)

        Nsize = 10
        arr = np.random.randint(0, 300, size=Nsize*2).reshape(Nsize, 2)
        d3 = fa.calc_distance(arr)
        d4 = traj.calc_distance(arr)
        d5 = pyca.calc_distance(traj, arr)
        d6 = pyca.calc_distance(fa, arr)
        d7 = pyca.calc_distance([fa, traj], arr, n_frames=2*fa.n_frames)
        aa_eq(d3, d4)
        aa_eq(d3, d5)
        aa_eq(d3, d6)
        aa_eq(d3, d7[:fa.n_frames])
        aa_eq(d3, d7[fa.n_frames:])


if __name__ == "__main__":
    unittest.main()
