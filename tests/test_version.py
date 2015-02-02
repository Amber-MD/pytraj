
import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        from pytraj.version import cppversion as cpptraj_version
        print("cpptraj_version = %s" % cpptraj_version) 

if __name__ == "__main__":
    unittest.main()
