from __future__ import print_function
import unittest
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists

# NOTE: no assert, just check for segfault


class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        trajiter = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        atm = traj.top("@CA")
        #print(traj[atm])
        #print(trajiter[atm])
        #print(traj[:5, atm])
        #print(trajiter[:5, atm])

        f0 = traj[5]
        #print(f0[atm])
        #print(traj[5, atm])
        #print(trajiter[5, atm])
        #print(traj[:, atm])
        #print(trajiter[:, atm])
        #print(traj[0, atm, 0])
        #print(trajiter[0, atm, 0])
        #print(traj[atm, 0, 0])
        #print(trajiter[atm, 0, 0])
        #print(traj[atm, 2:6, 0])
        #print(trajiter[atm, 2:6, 0])
        #print(traj[atm, :, 0])
        #print(trajiter[atm, :, 0])

        f0 = traj[0]
        f0.top = traj.top
        f0['@CA']
        traj[0, '@CA']


if __name__ == "__main__":
    unittest.main()
