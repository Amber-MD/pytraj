import unittest
from pytraj import *
import pytraj as pt
from utils import fn


class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))
        f0 = traj[0]
        # creat new frame based on f0 and atommask
        Frame(f0, traj.top('@CA'))

        atm = AtomMask("")
        traj.top._set_integer_mask(atm)


if __name__ == "__main__":
    unittest.main()
