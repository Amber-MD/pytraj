import unittest
import pytraj as pt
import sys
from pytraj.testing import aa_eq
from pytraj.core.c_core import Command
from pytraj.core.c_core import CpptrajState
from pytraj.compat import zip
from pytraj import iterframe_master
from pytraj import Frame, Trajectory

text = """
parm ./data/Tc5b.top
trajin ./data/Tc5b.x
rotate x 60 y 120 z 50 @CA
trajout rotated_frame0.x60y120z50.Tc5b.r
"""


def iter_me(obj, n_frames):
    it = iterframe_master(obj)
    for idx, frame in enumerate(it):
        pass
    assert idx + 1 == n_frames


class TestIterFrameMaster(unittest.TestCase):

    def test_iter(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        fa = traj[:]

        iter_me(traj, traj.n_frames)
        iter_me(fa, traj.n_frames)

        iter_me(traj(), traj.n_frames)
        iter_me(fa(), traj.n_frames)

        iter_me(traj(mask='@CA'), traj.n_frames)
        iter_me(fa(mask='@CA'), traj.n_frames)

        iter_me([traj, fa], 2 * traj.n_frames)
        iter_me((traj, fa), 2 * traj.n_frames)
        iter_me((traj, (fa[0], )), traj.n_frames + 1)

        for frame in iterframe_master(traj[0]):
            assert frame.n_atoms == traj.top.n_atoms

        i = 0
        for frame in iterframe_master([traj, traj[:1]]):
            i += 1
            assert frame.n_atoms == traj.top.n_atoms
        assert i == traj.n_frames + 1

        i = 0
        for frame in iterframe_master(traj.iterchunk()):
            i += 1
            assert isinstance(frame, Frame)
        assert i == traj.n_frames

        i = 0
        for frame in iterframe_master([traj.iterchunk(), ]):
            i += 1
            assert isinstance(frame, Frame)
        assert i == traj.n_frames

        # raise if wrong type
        def test_raise():
            for frame in pt.iterframe_master([0, 3]):
                pass

        self.assertRaises(TypeError, lambda: test_raise())

    def test_iter_with_a_list_of_frame_and_trajectory_and_FrameIterator(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        ref = traj[0]

        for idx, frame in enumerate(iterframe_master([traj[0], traj])):
            assert isinstance(frame, Frame), 'must a a Frame'
        assert idx == traj.n_frames

        # FrameIterator: traj()
        for idx, frame in enumerate(iterframe_master([traj[0], traj()])):
            assert isinstance(frame, Frame), 'must a a Frame'
        assert idx == traj.n_frames

        # FrameIterator: traj()
        fi = traj.iterframe(frame_indices=[0, 3])
        for idx, frame in enumerate(iterframe_master([traj[0], fi])):
            assert isinstance(frame, Frame), 'must a a Frame'

    def test_assert(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        fa = Trajectory.from_iterable(iterframe_master(traj), top=traj.top)

        for f0, f1 in zip(fa, traj):
            aa_eq(f0.xyz, f1.xyz)

    def test_TrajectorView(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        # make mutable traj
        t0 = traj[:]
        t1 = traj[:]
        indices = [2, 4, 6]
        # iter frame to return a view
        for f in t0.iterframe(frame_indices=indices):
            f.xyz += 1.0
        aa_eq(t0.xyz[indices], traj[indices].xyz + 1.)
        aa_eq(t1.xyz[indices], traj[indices].xyz)


class TestIterFrameFromArray(unittest.TestCase):

    def test_iterframe_from_array(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")

        # no mass
        fi = pt.iterframe_from_array(traj.xyz, traj.n_atoms,
                                     range(traj.n_frames))
        for f_traj, f_fi in zip(traj, fi):
            aa_eq(f_traj.xyz, f_fi.xyz)
            # f_fi mass must be [1.0, ...]
            aa_eq(f_fi.mass, [1.0 for _ in range(traj.n_atoms)])

        # mass from Topology
        fi = pt.iterframe_from_array(traj.xyz, traj.n_atoms,
                                     range(traj.n_frames), traj.top)
        for f_traj, f_fi in zip(traj, fi):
            aa_eq(f_traj.xyz, f_fi.xyz)
            aa_eq(f_fi.mass, f_traj.mass)
            aa_eq(f_fi.mass, traj.top.mass)

        # mass from array
        fi = pt.iterframe_from_array(traj.xyz,
                                     traj.n_atoms,
                                     range(traj.n_frames),
                                     mass=traj.top.mass)
        for f_traj, f_fi in zip(traj, fi):
            aa_eq(f_traj.xyz, f_fi.xyz)
            aa_eq(f_fi.mass, f_traj.mass)
            aa_eq(f_fi.mass, traj.top.mass)


if __name__ == "__main__":
    unittest.main()
