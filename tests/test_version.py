import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):
    def test_0(self):
        from pytraj.__cpptraj_version__ import __cpptraj_version__, info
        #print("cpptraj_version = %s" % __cpptraj_version__)
        #print(info())


if __name__ == "__main__":
    unittest.main()
