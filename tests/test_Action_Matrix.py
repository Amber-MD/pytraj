from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import allactions
from pytraj.cast_dataset import cast_dataset
from pytraj import adict

farray = TrajReadOnly(top=Topology("./data/Tc5b.top"), 
                    filename='data/md1_prod.Tc5b.x', 
                    )
class TestRadgyr(unittest.TestCase):
    def test_0(self):
        dslist = DataSetList()
        act = adict['matrix']()
        act.run(command="byres @CA", current_frame=farray, 
                current_top=farray.top, dslist=dslist)

        d1 = cast_dataset(dslist[0], dtype="matrix")
        print (d1.size)
        print (dir(d1))
        print (d1.n_cols, d1.n_rows)
        print (d1.dtype)
        print (d1.ndim)
        print (d1.kind)
        print (d1.data_format)
        # TODO : add assert to make sure reproducing cpptraj output

        for i in range(d1.size):
            print (d1[i])

if __name__ == "__main__":
    unittest.main()
