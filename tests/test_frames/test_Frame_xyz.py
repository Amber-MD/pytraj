from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):

    def test_0(self):
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        farray = traj[:]
        xyz_save = farray[0].xyz.copy()

        farray[9].xyz[:] = farray[0].xyz
        assert_almost_equal(farray[9].xyz, farray[0].xyz)
        assert_almost_equal(farray[9].xyz, xyz_save.flatten())

        f0 = Frame(traj.n_atoms)
        f0.xyz[:] = farray[0].xyz

        f1 = Frame()
        f1.append_xyz(farray[0].xyz)
        assert_almost_equal(f1.xyz, farray[0].xyz)


if __name__ == "__main__":
    unittest.main()
