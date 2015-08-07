import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj import adict


class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        f0 = traj[0]
        f0cp = f0.copy()
        act = adict['strip']
        act.read_input('!@CA', traj.top.copy())
        top = traj.top.copy()
        act.process(top)
        newframe = Frame()
        act.do_action(traj, newframe)


if __name__ == "__main__":
    unittest.main()
