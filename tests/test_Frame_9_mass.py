import unittest
import numpy as np
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        frame = traj[0].copy()
        print (frame.mass)

        print (frame.rmsd(traj[1], use_mass=False))

        frame.set_frame_m(traj.top)
        print (frame.mass)
        print (frame.rmsd(traj[1], use_mass=True))
        print (frame.rmsd(frame, use_mass=True))

        arr0 = []
        print (frame[0])
        for f0 in traj:
            print (f0[0])
            f0.set_frame_m(traj.top)
            arr0.append(frame.rmsd(f0, use_mass=True))

        # load cpptraj output
        rmsd_save = np.loadtxt("./data/rmsd_allatoms_to_1st.Tc5b.use_mass.dat", skiprows=1)
        rmsd_save = rmsd_save.transpose()
        print (rmsd_save[1])
        print (arr0)

        assert_almost_equal(arr0, rmsd_save[1])

    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        frame = traj[0].copy()
        print (frame.mass)
        frame.set_frame_m(traj.top)
        frame2 = Frame(frame, traj.top("@CA"))
        print (frame2.mass)

if __name__ == "__main__":
    unittest.main()
