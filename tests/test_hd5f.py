import unittest
import numpy as np
from pytraj.base import *
from pytraj import io as io
from pytraj.utils import aa_eq
from pytraj.decorators import test_if_having

class Test(unittest.TestCase):
    @test_if_having("h5py")
    def test_1(self):
        traj = io.load_hd5f("./data/ala2.h5")
        print (traj)
        assert traj.top.has_box() == False

    @test_if_having("h5py")
    @test_if_having("mdtraj")
    def test_1(self):
        import mdtraj as md
        from mdtraj.testing import get_fn
        fn = get_fn("frame0.h5")
        m_traj = md.load(fn)
        print (m_traj)
        traj = io.load_hd5f(fn)
        traj2 = io.load_hd5f(fn, autoconvert=False)
        print (traj)
        assert traj.top.has_box() == True
        assert traj2.top.has_box() == True

        aa_eq(traj.xyz, m_traj.xyz * 10)
        aa_eq(traj2.xyz, m_traj.xyz)
        print (traj.top.box, traj2.top.box)
        save_list10 = [10.0, 10.0, 10.0, 90.0, 90.0, 90.0]
        save_list = [1.0, 1.0, 1.0, 90.0, 90.0, 90.0]
        blist = traj.top.box.tolist()
        assert blist == save_list10

        for frame in traj:
            blist = frame.box.tolist()
            assert blist == save_list10

        # no autoconvert
        for frame in traj2:
            blist = frame.box.tolist()
            assert blist == save_list


if __name__ == "__main__":
    unittest.main()
