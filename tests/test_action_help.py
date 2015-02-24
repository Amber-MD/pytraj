import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        from pytraj import action_help as ah

        # print all keywords
        ah()

        # test specific keyword
        ah("rmsd")
        ah("spam")

if __name__ == "__main__":
    unittest.main()
