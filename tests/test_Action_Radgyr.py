from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import allactions
from pytraj.cast_dataset import cast_dataset
from pytraj import adict
from pytraj.common_actions import distance
import numpy as np
from rmsd import rmsd

farray = TrajReadOnly(top=Topology("./data/Tc5b.top"), 
                    filename='data/md1_prod.Tc5b.x', 
                    )

class TestRadgyr(unittest.TestCase):
    def test_0(self):
        dslist0 = DataSetList()
        dslist1 = DataSetList()

        act1 = adict['distance']()
        act1.run(":2@CA :10@CA", farray, farray.top, dslist1)
        act0 = adict['radgyr']()
        act0.run("@CA", farray, farray.top, dslist0)

        d0 = cast_dataset(dslist0[0])
        d1 = cast_dataset(dslist1[0])
        d0_0 = np.loadtxt("./data/radgyr.Tc5b.dat", skiprows=1).transpose()[1]
        d1_0 = np.loadtxt("./data/CAres2_CAres10.Tc5b.dat", skiprows=1).transpose()[1]

        print(d0.data[:10])
        print(d1.data[:10])
        print(d0_0[:10])
        print(d1_0[:10])

        # TODO : (fail)
        assert (rmsd(d1.data[:10], d1_0)) < 1E-3
        assert (rmsd(d0.data[:10], d0_0)) < 1E-3

        arr0 = np.empty(farray.size)
        for i, frame in enumerate(farray):
            frame.set_top(farray.top)
            #print ("mask of :2@CA", frame[':2@CA'])
            #print (type(frame[':2@CA']))
            arr0[i] = (distance(frame[':2@CA'][0], frame[':10@CA'][0]))

        assert (rmsd(d1_0, arr0)) < 1E-3

if __name__ == "__main__":
    unittest.main()
