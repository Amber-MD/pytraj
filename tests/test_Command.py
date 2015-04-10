import unittest
import sys
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.Command import Command
from pytraj.CpptrajState import CpptrajState

class Test(unittest.TestCase):
    def test_0(self):
        fname = "./tc5b.rotate.in"
        state = Command.get_state(fname)
        print (state)
        print (dir(state))
        print (state.toplist[0])

if __name__ == "__main__":
    unittest.main()
