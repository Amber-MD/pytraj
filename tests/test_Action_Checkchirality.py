import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj import allactions
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):
    def test_0(self):
        act = allactions.Action_CheckChirality()
        act.help()


if __name__ == "__main__":
    unittest.main()
