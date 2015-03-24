import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        f0 = traj[0]
        top = f0.get_top()
        if top:
            print("has top")
        else:
            print ("dont have top")
        print (top)
        f0.set_top(traj.top)
        top = f0.get_top()
        print (top)

        if top:
            print("has top")
        else:
            print ("dont have top")

if __name__ == "__main__":
    unittest.main()
