import unittest
import pytraj as pt
from pytraj import *
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal, is_generator
import numpy as np


class Test(unittest.TestCase):

    def test_0(self):
        dslist = DatasetList()
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")

        calc_distance = adict['distance']

        act = calc_distance
        act.read_input(":2@CA :10@CA", traj.top, dslist=dslist)
        act.setup(traj.top)
        act.compute((traj.iterchunk(2), ))
        assert dslist.size == 1
        assert dslist[0].size == traj.n_frames
        cppout = np.loadtxt("./data/CAres2_CAres10.Tc5b.dat",
                            skiprows=1).transpose()[1]
        assert_almost_equal(dslist[0][:], cppout)
        act.compute((traj.iterchunk(chunksize=4, stop=8), ))

    def test_1(self):
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        i = 0
        for farray in traj.iterchunk(chunksize=4, stop=8):
            i += 1
        assert i == 2

        i = 0
        for farray in traj.iterchunk(chunksize=4):
            i += 1
        assert i == 3

        i = 0
        for farray in traj.iterchunk(chunksize=2):
            i += 1
        assert i == 5

        i = 0
        for farray in traj.iterchunk(start=3, chunksize=4, stop=8):
            i += 1
        assert i == 2

        # action on chunk_iter

        pt.calc_distance(
            [traj.iterchunk(), traj.iterchunk(), traj[0]],
            '@CA @CB',
            top=traj.top)

        rmsd0 = pt.rmsd(traj.iterchunk(3), ref=traj[-1], top=traj.top)
        rmsd1 = pt.rmsd(traj, ref=-1)
        assert_almost_equal(rmsd0, rmsd1)


if __name__ == "__main__":
    unittest.main()
