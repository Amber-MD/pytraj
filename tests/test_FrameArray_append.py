import unittest
from pytraj._utils import _tease_FrameArray
from time import time
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import no_test
from pytraj import adict

class Test(unittest.TestCase):
    #@no_test
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        farray = FrameArray()

        t0 =  time()
        for i in range(1000):
            for frame in traj:
                farray.append(frame, copy=False)

        print (time() - t0)
        print (farray)
        farray.join(farray)
        print (farray)

    #@no_test
    def test_1(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        farray = FrameArray()
        farray.top = traj.top.copy()
        t0 = time()
        _tease_FrameArray(farray, traj[:], 10000)
        print (time()-t0)
        print (farray)

        farray2 = FrameArray()
        farray2.top = traj.top.copy()
        t0 = time()
        for frame in farray:
            farray2.append(frame, copy=False)
        print (time()-t0)

        t0 = time()
        for frame in farray:
            farray2.append(frame, copy=True)
        print (time()-t0)
        print (farray2)

        ## try doing action
        farray3 = farray2[:1000]
        farray3.join((farray3[:], farray3[:], farray3[:]))
        print (farray3)
        t0 = time()

        dslist = DataSetList()
        for f0 in farray3:
            assert f0.n_atoms == farray3.top.n_atoms
        print (farray3.top)
        #adict['distance'](":2@CA :10@CA", farray3,
        #                  farray3.top, dslist=dslist)
        
        d0 = adict['distance'](":2@CA :10@CA", farray3,
                               farray3.top, quick_get=True)
        print (time()-t0)
        print (d0.size)

if __name__ == "__main__":
    unittest.main()
