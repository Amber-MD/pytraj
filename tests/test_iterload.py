import unittest
from pytraj.io import iterload
from pytraj.base import *
from pytraj.externals.six import next

# test load single frame
# TODO : got Segmentation fault (core dumped) if run this script several times
# Reason: messed up with memoryview (FIXME) (see ./LOG/log_0.txt)
# Posible reason:
# firs, we assign a Frame view for farray[0] : frame0_0 = farray[0]
# then we assign farray to new value --> frame0 was deallocated because the old `farray` address was gone
# make copy?


class TestIterLoad(unittest.TestCase):
    def test_iterload_0(self):
        """get single frame"""
        genobj = iterload(filename="./data/md1_prod.Tc5b.x", top=Topology("./data/Tc5b.top"))
        
        frame0 = next(genobj)
        print(frame0)
        print(next(genobj))
        print(next(genobj)[0])
        print(next(genobj).n_atoms)

    def test_iterload_1(self):
        """get FrameArray"""
        print(__doc__)
        for i in range(1):
            print("")
            print("loop number = %s-th" % i)
            genobj = iterload(filename="./data/md1_prod.Tc5b.x", top=Topology("./data/Tc5b.top"), chunk=5)
            
            farray = next(genobj)
            farray.warning = True
            # create a view of farray[0]
            frame0_0 = farray[0].copy()
            print(farray)
            assert farray.size == 5
            print(farray[0].n_atoms)
            assert farray[0].n_atoms == 304

            print("next iteration")
            farray = next(genobj)
            frame0_1 = farray[0].copy()
            print(farray)
            assert farray.size == 5
            print(farray[0].n_atoms)
            assert farray[0].n_atoms == 304

            print("make sure that we DID the iteration")
            assert frame0_0.coords != frame0_1.coords
            #print frame0_0[:10]
            print(frame0_1.coords[:10])
            print(frame0_0.coords[:10])

    def test_iterload_2(self):
        # TODO : add more
        pass

if __name__ == "__main__":
    unittest.main()
