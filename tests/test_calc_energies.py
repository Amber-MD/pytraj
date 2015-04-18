from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca

class Test(unittest.TestCase):
    @test_if_having("sander")
    @test_if_having("chemistry")
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        edict = pyca.calc_energies("./data/Tc5b.top", traj)
        print (edict)

    @test_if_having("sander")
    @test_if_having("chemistry")
    @test_if_having("pandas")
    def test_1(self):
        from pytraj import to_dataframe
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        edict = pyca.calc_energies("./data/Tc5b.top", traj)
        print (to_dataframe(edict))

if __name__ == "__main__":
    unittest.main()
