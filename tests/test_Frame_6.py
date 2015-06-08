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

        # FIXME: the results are not matched yet. Need to check
        print (traj.top("@CA"))
        _rmsd_notfit = frame0.rmsd_nofit(ref, traj.top("@CA"))
        print ("_rmsd_notfit before fitting = ", _rmsd_notfit)
        frame0.rmsfit(ref, traj.top("@CA"))
        print (frame0.same_coords_as(frame0cp))
        print ('rmsd between frame0 and frame0cp after fit_to =', 
               frame0.rmsd(frame0cp))
        _rmsd_notfit = frame0.rmsd_nofit(ref, traj.top("@CA"))
        print ("_rmsd_notfit after using fit_to = ", _rmsd_notfit)
        _rmsd_fit = frame0.rmsd(ref, traj.top("@CA"))
        print ("_rmsd_fit after using fit_to = ", _rmsd_fit)

    def test_1(self):
        print ("test setitem for mask")
        import numpy as np
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        frame0 = traj[0]
        top = traj.top
        print (frame0[top('@CA')])
        atm = top("@CA")
        arr0 = np.arange(atm.n_atoms*3).reshape(atm.n_atoms, 3)
        frame0.set_top(top)
        print (frame0['@CA'][:] == arr0)
        print (frame0['@CA'][:])
        frame0[top('@CA')] = arr0
        print (frame0['@CA'][:])
        assert_almost_equal(frame0['@CA'][:].flatten(), arr0.flatten())

    def test_2(self):
        print ("test set Frame")
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        frame0 = traj[0]
        newframe = Frame(20)
        newframe.set_frame(frame0, traj.top('@CA'))
        frame0.strip_atoms('!@CA', traj.top)
        print (newframe[:10])
        print (frame0[:10])
        print (newframe.size)
        print (frame0.size)
        assert_almost_equal(frame0.coords, newframe.coords)

    def test_3(self):
        from pytraj.common_actions import calc_dihedral, calc_angle
        from pytraj.common_actions import calc_distance
        print ("test calc torsion, angle")
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        frame0 = traj[0]

        Nsize = 1000000
        np.random.seed(0)
        indices = np.random.randint(0, 300, size=Nsize*3).reshape(Nsize, 3)
        indices_dih = np.random.randint(0, 300, size=Nsize*4).reshape(Nsize, 4)
        indices_dist = np.random.randint(0, 300, size=Nsize*2).reshape(Nsize, 2)
        print (indices.shape)

        with Timer() as t:
            ang_0 = frame0.calc_angle(indices)
        print ("angle: time to calculate %s data points = %s (s)" % (Nsize, t.time_gap()))

        with Timer() as t:
            dih_0 = frame0.calc_dihedral(indices_dih)
        print ("dih: time to calculate %s data points = %s (s)" % (Nsize, t.time_gap()))

        with Timer() as t:
            dist_0 = frame0.calc_distance(indices_dist)
        print ("dist: time to calculate %s data points = %s (s)" % (Nsize, t.time_gap()))
        d0_saved = dist_0

        # randomly take 100 data points to assert
        for _ in range(10):
            idx = np.random.randint(0, Nsize-1)
            id0, id1 = indices_dist[idx] + 1
            dist_command = "@%s @%s" % (id0, id1)

            id0, id1, id2 = indices[idx] + 1
            angle_command = "@%s @%s @%s" % (id0, id1, id2)

            id0, id1, id2, id3 = indices_dih[idx] + 1
            dih_command = "@%s @%s @%s @%s" % (id0, id1, id2, id3)

            d0 = calc_dihedral(frame0, dih_command, traj.top, dtype='ndarray')
            d1 = calc_angle(frame0, angle_command, traj.top, dtype='ndarray')
            d3 = calc_distance(frame0, dist_command, traj.top, dtype='ndarray')
            print ('d0, d1, d3', d0, d1, d3)
            print (dist_0[idx], d3)
            print (ang_0[idx], d1)
            print (dih_0[idx], d0)


            SMALL = 1E-7
            from numpy import abs
            assert abs(dih_0[idx] - d0) < SMALL
            assert abs(ang_0[idx] - d1) < SMALL
            assert abs(dist_0[idx] - d3) < SMALL

if __name__ == "__main__":
    unittest.main()
