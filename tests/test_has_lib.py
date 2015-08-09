from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.utils import has_


class Test(unittest.TestCase):
    def test_0(self):
        if has_('pytraj'):
            print("has pytraj")

        if has_('numpy'):
            print("has numpy")


if __name__ == "__main__":
    unittest.main()
