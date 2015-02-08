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

if __name__ == "__main__":
    unittest.main()
