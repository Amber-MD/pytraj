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
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        trajin = """
        trajin ./data/md1_prod.Tc5b.x
        parm ./data/Tc5b.top
        reference ./data/Tc5b.crd
        trajout ./output/test_command.nc
        rms reference @CA out ./output/test_rmsd_command.dat
        """
        state = CpptrajState()
        print (dir(Command))
        Command.process_input(state, trajin)
        print (dir(state))
        print (state.is_empty())

if __name__ == "__main__":
    unittest.main()
