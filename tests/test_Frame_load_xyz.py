from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import assert_almost_equal, Timer
from pytraj.decorators import test_if_having


class Test(unittest.TestCase):
    @Timer()
    @test_if_having("numpy")
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        frame = Frame()
        frame.append_xyz(traj[0].xyz)
        print(frame.coords[:10])
        print(traj[0].coords[:10])
        assert_almost_equal(frame.coords, traj[0].coords)
        assert_almost_equal(frame.xyz.flatten(), traj[0].xyz.flatten())
        print(frame.buffer2d.shape)
        print(traj[0].xyz.shape)
        assert_almost_equal(frame.buffer1d, traj[0].xyz.flatten())

    @Timer()
    @test_if_having("numpy")
    @test_if_having("mdtraj")
    def test_1(self):
        from mdtraj.formats import psf
        import mdtraj as md
        from mdtraj.testing import get_fn
        import numpy as np

        fname = get_fn('ala_ala_ala.pdb')
        m_traj = md.load(fname)
        print(m_traj)
        f0 = Frame()
        f1 = f0.copy()
        f0.append_xyz(m_traj.xyz[0].astype(np.float64))
        farray = mdio.load_mdtraj(m_traj, autoconvert=False, top=fname)
        f1 = farray[0]

        assert_almost_equal(f0.coords, f1.coords)


if __name__ == "__main__":
    unittest.main()
