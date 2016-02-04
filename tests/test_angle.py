from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.testing import cpptraj_test_dir


class Test(unittest.TestCase):

    def test_0(self):
        import numpy as np
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        fa = traj[:]
        mask = ':1@CA :14@CB :15CA'
        d0 = pt.calc_angle(traj, mask, dtype='dataset').to_ndarray()
        d1 = pt.angle(traj, mask)
        d2 = pt.angle(fa, mask)

        aa_eq(d0, d1)
        aa_eq(d0, d2)

        Nsize = 10
        arr = np.random.randint(0, 300, size=Nsize * 3).reshape(Nsize, 3)
        d3 = pt.angle(fa, arr)
        d4 = pt.angle(traj, arr)
        d5 = pt.calc_angle(traj, arr)
        d6 = pt.calc_angle(fa, arr)
        d7 = pt.calc_angle([fa, traj], arr, n_frames=2 * fa.n_frames)
        aa_eq(d3, d4)
        aa_eq(d3, d5)
        aa_eq(d3, d6)
        aa_eq(d3.T, d7.T[:fa.n_frames])
        aa_eq(d3.T, d7.T[fa.n_frames:])

        d8 = pt.angle(traj, mask, dtype='dataset')
        d9 = pt.tools.dict_to_ndarray(pt.angle(traj, mask, dtype='dict'))
        aa_eq(d0, d8.values)
        aa_eq(d0, d9)


if __name__ == "__main__":
    unittest.main()
