import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        from pytraj.actions.Action_Mask_2 import Action_Mask_2
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        act = Action_Mask_2()
        act.help()

if __name__ == "__main__":
    unittest.main()
