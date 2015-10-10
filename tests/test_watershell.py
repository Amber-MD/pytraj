from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.action_dict import ADICT
from pytraj.actions.CpptrajActions import Action_Watershell


class TestWatershell(unittest.TestCase):
    def test_watershell(self):
        traj = pt.iterload("data/tz2.truncoct.nc",
                           "data/tz2.truncoct.parm7")
        d0 = pt.watershell(traj, '!:WAT')

if __name__ == "__main__":
    unittest.main()
