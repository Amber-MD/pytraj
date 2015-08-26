import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import _import_numpy


class TestHasnumpy(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        frame0 = traj[0]
        has_numpy, _np = _import_numpy()
        #print(has_numpy)
        #print(_np)


if __name__ == "__main__":
    unittest.main()
