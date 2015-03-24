from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    # TODO : add assert
    # what cpptraj dump in dslist?
    def test_0(self):
        fn = "./data/Test_NAstruct/adh026.3.pdb"
        traj = mdio.load(fn, fn)
        print (traj)

        dslist = DataSetList()
        act = adict['nastruct']
        act("", traj, dslist=dslist)
        act.print_output()
        print (dslist.size)
        
        d0 = dslist[0]
        for d0 in dslist:
            print (d0)

        dsmall = dslist.get_dataset(dtype='float')
        ds_int = dslist.get_dataset(dtype='integer')

        import numpy as np
        dsmall = np.asarray(dsmall)
        print (dsmall.shape)

        print (dslist.get_legends())
        print (dslist['2G15C-3G14C'][0][:])
        print (dslist['1G16C'])
        print (dslist['1G16C'].__len__())
        print ('traj.n_frames = %s' % traj.n_frames)

        print (ds_int)

if __name__ == "__main__":
    unittest.main()
