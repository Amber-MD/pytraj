import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        print (traj['@CA'])
        frame0 = traj[0]
        print (hasattr(frame0, 'shape'))
        assert_almost_equal(frame0[traj.top("@CA")].flatten(), 
                            traj['@CA'][0].flatten())

if __name__ == "__main__":
    unittest.main()
