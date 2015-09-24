from __future__ import print_function
import pytraj as pt
import unittest
import pytraj as pt
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.testing import cpptraj_test_dir, local_test
import pytraj.common_actions as pyca
from pytraj.misc import from_legends_to_indices
from pytraj.utils import Timer
import numpy as np


class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        fa = traj[:]
        mask = ':1@CA :14@CB'
        d0 = pyca.calc_distance(traj, mask)
        d1 = pt.distance(traj, mask)
        d2 = pt.calc_distance(fa, mask)

        aa_eq(d0, d1)
        aa_eq(d0, d2)

        Nsize = 12
        arr = np.random.randint(0, 300, size=Nsize * 2).reshape(Nsize, 2)
        d3 = pt.calc_distance(fa, arr)
        d4 = pt.distance(traj, arr)
        d5 = pyca.calc_distance(traj, arr)
        d6 = pyca.calc_distance(fa, arr)
        d7 = pyca.calc_distance([fa, traj], arr, n_frames=2 * fa.n_frames)
        aa_eq(d3, d4)
        aa_eq(d3, d5)
        aa_eq(d3, d6)
        aa_eq(d3.T, d7.T[:fa.n_frames])
        aa_eq(d3.T, d7.T[fa.n_frames:])

    def test_2(self):
        # calculate distance without specifying n_frames
        # TrajectoryIterator
        import numpy as np
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        arr = pt.distance(traj(stop=4), [0, 5])
        arr1 = pt.distance(traj(stop=4), [0, 5], n_frames=4)
        #print(arr, arr1)
        assert np.all(arr == arr1)

        arr2 = pt.distance(traj(stop=1000), [0, 5])
        arr3 = pt.distance(traj(stop=traj.n_frames), [0, 5])
        assert np.all(arr2 == arr3)

        # Trajectory
        traj = traj[:]
        arr = pt.distance(traj(stop=4), [0, 5])
        arr1 = pt.distance(traj(stop=4), [0, 5], n_frames=4)
        assert np.all(arr == arr1)

        arr2 = pt.distance(traj(stop=1000), [0, 5])
        arr3 = pt.distance(traj(stop=traj.n_frames), [0, 5])
        assert np.all(arr2 == arr3)


if __name__ == "__main__":
    unittest.main()
