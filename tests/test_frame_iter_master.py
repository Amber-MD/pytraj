import unittest
import sys
from pytraj.base import *
from pytraj import Frame
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.Command import Command
from pytraj.CpptrajState import CpptrajState
from pytraj.decorators import no_test
from pytraj.six_2 import izip as zip

text="""
parm ./data/Tc5b.top
trajin ./data/md1_prod.Tc5b.x
rotate x 60 y 120 z 50 @CA
trajout rotated_frame0.x60y120z50.Tc5b.r
"""

def iter_me(obj, n_frames):
    from pytraj._shared_methods import _frame_iter_master
    it = _frame_iter_master(obj) 
    for idx, frame in enumerate(it):
        pass
    assert idx + 1 == n_frames

class Test(unittest.TestCase):
    # Potentail failing: using iterator of iterator
    # for example: it = _frame_iter_master(traj)
    # for idx, frame in enumerate(it):
    #     traj[idx]
    # will give segmentation fault
    @no_test
    def test_cpptraj_file(self):
        # FIXME: got segmentation fault.
        # reason: when iterating TrajinLis, no information about Topology in `traj`
        # --> can not read frame correctly
        from pytraj._shared_methods import _frame_iter_master
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        fname = "./tc5b.rotate.in"
        with open(fname, 'w') as f:
            f.write(text)
        state = mdio.load_cpptraj_file(fname)
        print (state)
        trajinlist = state.get_trajinlist()

        for idx, frame in enumerate(_frame_iter_master(trajinlist)):
            pass
        assert idx + 1 == traj.n_frames

    def test_0(self):
        from pytraj._shared_methods import _frame_iter_master
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        fa = traj[:]

        print ("iter traj")
        iter_me(traj, traj.n_frames)
        iter_me(fa, traj.n_frames)

        print ("iter of traj frame_iter")
        iter_me(traj(), traj.n_frames)
        iter_me(fa(), traj.n_frames)

        print ("iter of traj frame_iter with mask")
        iter_me(traj(mask='@CA'), traj.n_frames)
        iter_me(fa(mask='@CA'), traj.n_frames)

        print ("iter list/tuple")
        iter_me([traj, fa], 2 * traj.n_frames)
        iter_me((traj, fa), 2 * traj.n_frames)
        iter_me((traj, (fa[0],)), traj.n_frames + 1)

        print ("iter frame")
        for frame in _frame_iter_master(traj[0]):
            assert frame.n_atoms == traj.top.n_atoms

        print ("iter frame")
        i = 0
        for frame in _frame_iter_master([traj, traj[:1]]):
            i += 1
            assert frame.n_atoms == traj.top.n_atoms
        assert i == traj.n_frames + 1

        print ("iter chunk_iter")
        i = 0
        for frame in _frame_iter_master(traj.chunk_iter()):
            i += 1
            assert isinstance(frame, Frame)
        assert i == traj.n_frames

        print ("list of chunk_iter")
        i = 0
        for frame in _frame_iter_master([traj.chunk_iter(),]):
            i += 1
            assert isinstance(frame, Frame)
        assert i == traj.n_frames

    def test_assert(self):
        from pytraj._shared_methods import _frame_iter_master as _it_f
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        fa = Trajectory()
        fa.top = traj.top.copy()
        fa.load(_it_f(traj))

        for f0, f1 in zip(fa, traj):
            print (f0[0, :], f1[0, :])
            assert_almost_equal(f0.coords, f1.coords)

if __name__ == "__main__":
    unittest.main()
