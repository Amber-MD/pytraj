import unittest

from pytraj import io as mdio


class Test(unittest.TestCase):

    def test_0(self):
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")[:]
        arr0 = traj[:, :, :]
        arr0[0, 0, 0] = 105.


if __name__ == "__main__":
    unittest.main()
