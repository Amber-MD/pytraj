import unittest
from pytraj.base import *
from pytraj.utils.check_and_assert import assert_almost_equal, eq

class TestAtomMask(unittest.TestCase):
    def test_0(self):
        atm = AtomMask("@CA")
        print(atm.n_atoms)
        top = Topology("./data/Tc5b.top")
        print(dir(atm))
        print(dir(top))
        top.set_integer_mask(atm)
        print(atm.n_atoms)
        atm.brief_mask_info()
        print(atm.mask_string)

    def test_1(self):
        print("test_1")
        atm = AtomMask(10)
        atm.brief_mask_info()
        print("end test_1")

    def test_2(self):
        print("heavy")
        print("Oxygen")
        print("hydrogen")

    def test_3_indexing(self):
        print("test_indexing")
        top = Topology("./data/Tc5b.top")
        atm = AtomMask("@CA")
        top.set_integer_mask(atm)
        print(atm[0])

        for i in atm:
            print(i)

        print(dir(atm))

    def test_4(self):
        from array import array
        print ("add array")
        atm = AtomMask()
        indices = array('i', range(100))
        atm.add_selected_indices(indices)
        assert_almost_equal(indices, atm.selected_indices())

    def test_5_speed(self):
        top = Topology("./data/DOPC.parm7")
        from time import time

        t0 = time()
        indices = top(":WAT").indices
        gap_0 = time() - t0
        print ("time to call indices = %s (s)" % gap_0)

        t0 = time()
        _indices_view = top(":WAT")._indices_view
        gap_1 = time() - t0
        print ("time to call indices = %s (s)" % gap_1)

        print ("speed up when using memview = %s" % (gap_0/gap_1))
        count = 0 
        for i, j in zip(indices, _indices_view):
            if not i == j:
                count += 1
                print (i, j)
        assert count == 0

    def test_6_speed(self):
        import numpy as np
        from pytraj import AtomMask
        # test constructor from list/array/python array
        top = Topology("./data/DOPC.parm7")
        indices = top.select(":WAT")

        atm1 = AtomMask(indices)
        atm2 = AtomMask(list(indices))
        atm3 = AtomMask(np.asarray(indices))

        assert atm1.indices == atm2.indices == atm3.indices
        
if __name__ == "__main__":
    unittest.main()

