import unittest
import sys
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.Command import Command
from pytraj.CpptrajState import CpptrajState

parm  = "./data/tz2.truncoct.parm7"
trajin = "./data/tz2.truncoct.nc"
gist_command = "doorder doeij gridcntr 25.0 31.0 30.0 griddim 41 41 45 gridspacn 0.50"

class Test(unittest.TestCase):

    def test_1(self):
        from pytraj import calculate
        traj = mdio.load(trajin, parm)
        print (traj)
        dslist = calculate("gist", gist_command, traj[0], top=traj.top)
        print (dslist.get_legends())

if __name__ == "__main__":
    unittest.main()
