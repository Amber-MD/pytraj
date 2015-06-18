from __future__ import print_function
import unittest
from pytraj import io
from pytraj.utils import eq, aa_eq

class Test(unittest.TestCase):
    def test_0(self):
        from pytraj import TrajectoryIterator
        # normal constructor
        traj = TrajectoryIterator("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        print (traj)
        n_frame0 = 10
        assert traj.n_frames == n_frame0

        # load a list of files
        N = 3
        traj = TrajectoryIterator(["./data/md1_prod.Tc5b.x" for _ in range(N)],
                                   "./data/Tc5b.top")
        print (traj)
        assert traj.n_frames == n_frame0 * N

        traj = TrajectoryIterator("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        print (traj)
        n_frame0 = 10
        assert traj.n_frames == n_frame0

        # a list of files with frame_slice argument
        traj = TrajectoryIterator(["./data/md1_prod.Tc5b.x" for _ in range(N)],
                                   "./data/Tc5b.top", frame_slice=[(0, 5), (3, 8)])
        print (traj)
        assert traj.n_frames == 20

        # a list of files with frame_slice argument but not using `frame_slice = `
        traj = TrajectoryIterator(["./data/md1_prod.Tc5b.x" for _ in range(N)],
                                   "./data/Tc5b.top", [(0, 5), (3, 8)])
        print (traj)
        assert traj.n_frames == 20


        # dummy constructor without Topology. need to catch ValueError
        self.assertRaises(ValueError, lambda : TrajectoryIterator(["./data/md1_prod.Tc5b.x" 
                                               for _ in range(N)]))

        # constructor with a Topology
        traj = TrajectoryIterator(top="./data/Tc5b.top")
        assert not traj.top.is_empty()
        assert traj.top.n_atoms == 304
        print (traj.top)

        # dummy constructor with one arg but don't know what it is
        # raise ValueError
        self.assertRaises(ValueError, lambda : TrajectoryIterator("./data/Tc5b.top"))

        # dummy constructor with frame_slice
        # TODO: asserWarns?
        #self.assertWarns(UserWarning, TrajectoryIterator(frame_slice=(0, 4)))
        TrajectoryIterator(frame_slice=(0, 4))

    def test_1(self):
        # test load using `pytraj.io`
        import pytraj as pt

        traj = pt.io.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        print (traj)
        n_frame0 = 10
        assert traj.n_frames == n_frame0

        # load a list of files
        N = 3
        traj = pt.io.iterload(["./data/md1_prod.Tc5b.x" for _ in range(N)],
                                   "./data/Tc5b.top")
        print (traj)
        assert traj.n_frames == n_frame0 * N

        traj = pt.io.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        print (traj)
        n_frame0 = 10
        assert traj.n_frames == n_frame0

        # a list of files with frame_slice argument
        traj = pt.io.iterload(["./data/md1_prod.Tc5b.x" for _ in range(N)],
                                   "./data/Tc5b.top", frame_slice=[(0, 5), (3, 8)])
        print (traj)
        assert traj.n_frames == 20

        # a list of files with frame_slice argument but not using `frame_slice = `
        self.assertRaises(TypeError, lambda : pt.io.iterload(["./data/md1_prod.Tc5b.x" for _ in range(N)],
                               "./data/Tc5b.top", [(0, 5), (3, 8)]))

if __name__ == "__main__":
    unittest.main()
