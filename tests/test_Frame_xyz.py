from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import no_test, test_if_having

class Test(unittest.TestCase):
    @test_if_having("numpy")
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        farray = traj[:]
        xyz_save = farray[0].xyz.copy()

        farray[9].xyz[:] = farray[0].xyz
        assert_almost_equal(farray[9].coords, farray[0].coords)
        assert_almost_equal(farray[9].coords, xyz_save.flatten())

        f0 = Frame(traj.n_atoms)
        f0.xyz[:] = farray[0].xyz

        f1 = Frame()
        f1.append_xyz(farray[0].xyz)
        assert_almost_equal(f1.coords, farray[0].coords)

if __name__ == "__main__":
    unittest.main()
