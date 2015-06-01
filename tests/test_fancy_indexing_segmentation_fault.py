from __future__ import print_function
import unittest
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists

# NOTE: no assert, just check for segfault
class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        atm = traj.top("@CA")
        print (traj[atm])
        print (traj[:5, atm])

        f0 = traj[5]
        print (f0[atm])
        print (traj[5, atm])
        print (traj[:, atm])
        print (traj[0, atm, 0])
        print (traj[atm, :, 0])

if __name__ == "__main__":
    unittest.main()
