import unittest
from pytraj.base import *
from pytraj.utils import assert_almost_equal

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

    def test_indexing(self):
        print("test_indexing")
        top = Topology("./data/Tc5b.top")
        atm = AtomMask("@CA")
        top.set_integer_mask(atm)
        print(atm[0])

        for i in atm:
            print(i)

        print(dir(atm))

    def test_3(self):
        from array import array
        print ("add array")
        atm = AtomMask()
        indices = array('i', range(100))
        atm.add_selected_indices(indices)
        assert_almost_equal(indices, atm.selected_indices())
        

if __name__ == "__main__":
    unittest.main()

