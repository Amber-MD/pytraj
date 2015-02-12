import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        frame0 = traj[0]
        frame0cp = frame0.copy()
        ref = traj[1]

        # FIXME: the results are not matched yet. Need to check
        print (traj.top("@CA"))
        _rmsd_notfit = frame0.rmsd_nofit(ref, traj.top("@CA"))
        print ("_rmsd_notfit before fitting = ", _rmsd_notfit)
        frame0.fit_to(ref, traj.top("@CA"))
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
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        frame0 = traj[0]
        top = traj.top
        print (frame0[top('@CA')])
        atm = top("@CA")
        arr0 = np.arange(atm.n_selected*3).reshape(atm.n_selected, 3)
        frame0.set_top(top)
        print (frame0['@CA'][:] == arr0)
        print (frame0['@CA'][:])
        frame0[top('@CA')] = arr0
        print (frame0['@CA'][:])
        assert_almost_equal(frame0['@CA'][:].flatten(), arr0.flatten())

    def test_2(self):
        print ("test set Frame")
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
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
        print ("test calc torsion, angle")
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        frame0 = traj[0]
        dih_0 = frame0.calc_dihedral(1, 2, 3, 4)
        ang_0 = frame0.calc_angle(1, 2, 3)
        d0 = calc_dihedral("@2 @3 @4 @5", frame0, traj.top)
        d1 = calc_angle("@2 @3 @4 @5", frame0, traj.top)
        assert dih_0 == d0[0]
        assert ang_0 == d1[0]

if __name__ == "__main__":
    unittest.main()
