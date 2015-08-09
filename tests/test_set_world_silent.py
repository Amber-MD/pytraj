from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):
    def test_0(self):
        from pytraj import set_world_silent

        print("no verbose")
        set_world_silent(True)
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        print("too many words")
        set_world_silent(False)
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        #
        traj.top.bond_info()


if __name__ == "__main__":
    unittest.main()
