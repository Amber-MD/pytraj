import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.action_dict import ActionDict

class Test(unittest.TestCase):
    def test_0(self):
        adict = ActionDict()
        print (adict)
        print (adict['rmsd'])

if __name__ == "__main__":
    unittest.main()
