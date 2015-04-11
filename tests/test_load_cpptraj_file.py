import unittest
import sys
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.Command import Command
from pytraj.CpptrajState import CpptrajState

text="""
parm ./data/Tc5b.top
trajin ./data/md1_prod.Tc5b.x
rotate x 60 y 120 z 50 @CA
trajout rotated_frame0.x60y120z50.Tc5b.r
"""

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        fname = "./tc5b.rotate.in"
        with open(fname, 'w') as f:
            f.write(text)
        state = Command.get_state(fname)
        assert (state.is_empty() == False)
        print (state)
        print (dir(state))
        assert (state.toplist[0].n_atoms == traj.top.n_atoms)
        _ctraj = state.get_trajinlist()[0]
        _ctraj.top = traj.top.copy()
        assert _ctraj.n_frames == traj.n_frames
        from pytraj import set_world_silent

        # turn off set_world_silent to see the cpptraj's talk
        set_world_silent(False)
        state.run()

    def test_1(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        fname = "./tc5b.rotate.in"
        with open(fname, 'w') as f:
            f.write(text)
        state = mdio.load_cpptraj_file(fname)
        assert (state.is_empty() == False)
        print (state)
        print (dir(state))
        assert (state.toplist[0].n_atoms == traj.top.n_atoms)
        _ctraj = state.get_trajinlist()[0]
        _ctraj.top = traj.top.copy()
        assert _ctraj.n_frames == traj.n_frames

if __name__ == "__main__":
    unittest.main()
