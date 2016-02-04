from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.utils import Timer


class Test(unittest.TestCase):

    def test_0(self):
        with Timer() as t:
            traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")


if __name__ == "__main__":
    unittest.main()
