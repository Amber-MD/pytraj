import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        top = mdio.loadpdb_rcsb("1l1h")
        print(top)

if __name__ == "__main__":
    unittest.main()
