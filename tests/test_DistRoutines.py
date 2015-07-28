import math
import unittest
from pytraj.base import *
from pytraj.math import Vec3
from pytraj.math.DistRoutines import distance
from pytraj.decorators import no_test
from pytraj import allactions


class Test(unittest.TestCase):

    def test_Vec3_noimage(self):
        arr0 = [0., 0., 1.]
        arr1 = [0., 0., 1.]
        v1 = Vec3(0., 0., 1.)
        v2 = Vec3(0., 0., 1.)
        v3 = Vec3(0., 0., 2.)
        v4 = Vec3(1., 1., 2.)
        assert distance(v1, v2) == 0
        assert distance(v1, v3) == 1.
        assert distance(v1, v4) == math.sqrt(1 + 1 + 1)

    def test_list_or_arr_noimage(self):
        arr0 = [0., 0., 1.]
        arr1 = [0., 0., 1.]
        arr4 = (1., 1., 2.)
        assert distance(arr0, arr1) == 0
        assert distance(arr0, arr4) == math.sqrt(1 + 1 + 1)

    def test_distance_from_frame(self):
        traj = Trajectory(
            filename="./data/md1_prod.Tc5b.x", top="./data/Tc5b.top")
        frame0 = traj[9]

        # print distance(frame0.coords[0:3], frame0.coords[96:99])
        for i in range(frame0.n_atoms):
            for j in range(i, frame0.n_atoms):
                distance(frame0.atoms(i), frame0.atoms(j))

if __name__ == "__main__":
    unittest.main()
