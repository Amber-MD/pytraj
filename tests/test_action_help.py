import unittest
from pytraj import info
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):

    def test_0(self):
        from pytraj import info
        info(adict["rmsd"])
        info(adict["spam"])

if __name__ == "__main__":
    unittest.main()
