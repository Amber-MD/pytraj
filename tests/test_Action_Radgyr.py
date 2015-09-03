from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import allactions
from pytraj.datasets import cast_dataset
from pytraj import adict
from pytraj.common_actions import distance
import numpy as np
from rmsd import rmsd

farray = TrajectoryIterator(top=Topology("./data/Tc5b.top"),
                            filename='data/md1_prod.Tc5b.x', )


class TestRadgyr(unittest.TestCase):
    def test_0(self):
        dslist0 = DataSetList()
        dslist1 = DataSetList()

        act1 = adict['distance']
        act1(":2@CA :10@CA", farray, farray.top, dslist1)
        act0 = adict['radgyr']
        act0("@CA", farray, farray.top, dslist0)

        d0 = cast_dataset(dslist0[0])
        d1 = cast_dataset(dslist1[0])
        d0_0 = np.loadtxt("./data/radgyr.Tc5b.dat", skiprows=1).transpose()[1]
        d1_0 = np.loadtxt(
            "./data/CAres2_CAres10.Tc5b.dat",
            skiprows=1).transpose()[1]


        assert (rmsd(d1.data[:10], d1_0)) < 1E-3
        assert (rmsd(d0.data[:10], d0_0)) < 1E-3

        arr0 = np.empty(farray.n_frames)
        for i, frame in enumerate(farray):
            frame.top = farray.top
            arr0[i] = distance(frame[':2@CA'][0], frame[':10@CA'][0])

        assert (rmsd(d1_0, arr0)) < 1E-3


if __name__ == "__main__":
    unittest.main()
