import unittest
from pytraj.base import *

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
        print(atm.mask_string())

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

if __name__ == "__main__":
    unittest.main()

