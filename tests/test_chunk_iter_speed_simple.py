from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq, eq_coords
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir, duplicate_traj
import pytraj.common_actions as pyca

class Test(unittest.TestCase):
    def test_1(self):
        from pytraj.testing import make_fake_traj
        from pytraj.compat import zip

        fa = make_fake_traj(1000, 1000)
        fa2 = fa.__class__(top=fa.top)

        for chunk in fa.chunk_iter(chunksize=100):
            fa2.join(chunk, copy=False)

        print (fa2.n_frames, fa.n_frames)
        assert fa2.n_frames == fa.n_frames
        eq_coords(fa2, fa)


        from pytraj.utils import Timer
        @Timer()
        def test_chunk(chunksize=2): 
            print ("chunk = %s" % chunksize)
            for chunk in fa.chunk_iter(chunksize=chunksize):
                pass

        @Timer()
        def test_frame_iter():
            for frame in fa:
                pass

        print ("test_chunk")
        test_chunk(2)
        test_chunk(10)
        test_chunk(100)
        test_chunk(200)
        test_chunk(500)
        print ("test frame_iter")
        test_frame_iter()

    def test_0(self):
        fa = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        f_list = []

        for chunk in fa.chunk_iter(chunksize=2):
            f_list.append(chunk.average())

        print (f_list)
        f0 = fa[:2].average()
        aa_eq(f_list[0].xyz, f0.xyz)

        f0 = fa[2:4].average()
        aa_eq(f_list[1].xyz, f0.xyz)

if __name__ == "__main__":
    unittest.main()
