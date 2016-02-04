import unittest
import numpy as np
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):

    def test_0(self):
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        frame = traj[0].copy()

        frame.set_mass(traj.top)

        arr0 = []
        for f0 in traj:
            f0.set_mass(traj.top)
            arr0.append(frame.rmsd(f0, use_mass=True))

        # load cpptraj output
        rmsd_save = np.loadtxt("./data/rmsd_allatoms_to_1st.Tc5b.use_mass.dat",
                               skiprows=1)
        rmsd_save = rmsd_save.transpose()

        assert_almost_equal(arr0, rmsd_save[1])

    def test_0(self):
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        frame = traj[0].copy()
        frame.set_mass(traj.top)
        frame2 = Frame(frame, traj.top("@CA"))


if __name__ == "__main__":
    unittest.main()
