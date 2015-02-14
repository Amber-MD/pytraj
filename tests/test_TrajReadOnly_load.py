import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        traj = TrajReadOnly("./data/md1_prod.Tc5b.x", "./data/Tc5b.top", arglist=ArgList('1 8 2'))
        print (traj)

if __name__ == "__main__":
    unittest.main()
