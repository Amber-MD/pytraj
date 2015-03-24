from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        from pytraj import set_world_silent
        set_world_silent(False)
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        f0 = traj[0]
        act = adict['makestructure']
        act("alpha:2-19", f0, traj.top)
        f0.save("./output/test_make_structure.pdb", traj.top, overwrite=True)

if __name__ == "__main__":
    unittest.main()
