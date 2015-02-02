import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.misc import action_dict

class TestMisc(unittest.TestCase):
    def test_0(self):
        print(list(action_dict.keys()))

if __name__ == "__main__":
    unittest.main()
