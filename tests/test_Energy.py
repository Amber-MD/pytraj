import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.Energy import Energy_Amber

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        ene = Energy_Amber()
        print (ene.E_bond(traj[0], traj.top, traj.top('*')))
        print (ene.E_angle(traj[0], traj.top, traj.top('*')))

if __name__ == "__main__":
    unittest.main()
