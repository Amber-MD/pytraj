import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
import numpy as np
from pytraj.utils import Timer


class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        frame0 = traj[0]
        frame0cp = frame0.copy()
        ref = traj[1]

        _rmsd_notfit = frame0.rmsd_nofit(ref, traj.top("@CA"))
        frame0.rmsfit(ref, traj.top("@CA"))
        _rmsd_notfit = frame0.rmsd_nofit(ref, traj.top("@CA"))
        _rmsd_fit = frame0.rmsd(ref, traj.top("@CA"))

    def test_1(self):
        import numpy as np
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        frame0 = traj[0]
        top = traj.top
        atm = top("@CA")
        arr0 = np.arange(atm.n_atoms * 3).reshape(atm.n_atoms, 3)
        frame0.top = top
        frame0[top('@CA')] = arr0
        assert_almost_equal(frame0['@CA'][:].flatten(), arr0.flatten())

    def test_2(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        frame0 = traj[0]
        newframe = Frame(20)
        newframe.set_frame(frame0, traj.top('@CA'))
        frame0.strip_atoms('!@CA', traj.top)
        assert_almost_equal(frame0.coords, newframe.coords)

    def test_3(self):
        from pytraj.common_actions import calc_dihedral, calc_angle
        from pytraj.common_actions import calc_distance
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        frame0 = traj[0]

        Nsize = 1000000
        np.random.seed(0)
        indices = np.random.randint(0, 300, size=Nsize * 3).reshape(Nsize, 3)
        indices_dih = np.random.randint(
            0, 300,
            size=Nsize * 4).reshape(Nsize, 4)
        indices_dist = np.random.randint(
            0, 300,
            size=Nsize * 2).reshape(Nsize, 2)

        with Timer() as t:
            ang_0 = frame0.calc_angle(indices)

        with Timer() as t:
            dih_0 = frame0.calc_dihedral(indices_dih)

        with Timer() as t:
            dist_0 = frame0.calc_distance(indices_dist)
        d0_saved = dist_0

        # randomly take 100 data points to assert
        for _ in range(10):
            idx = np.random.randint(0, Nsize - 1)
            id0, id1 = indices_dist[idx] + 1
            dist_command = "@%s @%s" % (id0, id1)

            id0, id1, id2 = indices[idx] + 1
            angle_command = "@%s @%s @%s" % (id0, id1, id2)

            id0, id1, id2, id3 = indices_dih[idx] + 1
            dih_command = "@%s @%s @%s @%s" % (id0, id1, id2, id3)

            d0 = calc_dihedral(frame0, dih_command, traj.top, dtype='ndarray')
            d1 = calc_angle(frame0, angle_command, traj.top, dtype='ndarray')
            d3 = calc_distance(frame0, dist_command, traj.top, dtype='ndarray')

            SMALL = 1E-7
            from numpy import abs
            assert abs(dih_0[idx] - d0) < SMALL
            assert abs(ang_0[idx] - d1) < SMALL
            assert abs(dist_0[idx] - d3) < SMALL


if __name__ == "__main__":
    unittest.main()
