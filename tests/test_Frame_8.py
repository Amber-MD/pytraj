import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):
    def test_0(self):
        print("test make Frame")
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        f0 = traj[0]
        # creat new frame based on f0 and atommask
        f1 = Frame(f0, traj.top('@CA'))
        print(f1)

        atm = AtomMask("")
        traj.top.set_integer_mask(atm)
        print(atm.n_atoms)


if __name__ == "__main__":
    unittest.main()
