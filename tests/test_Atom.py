import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        at = Atom()
        
        print(at)
        print(at.atomic_number)
        
        top = Topology("./data/HP36.top")
        a0 = top.atomlist[0]
        print(a0)
        print(a0.element)
        print(len(a0.element_short_name))
        print(a0.atomic_number)
        print(dir(a0))
        print(a0.molnum)
        print(a0.nbonds)
        print(a0.typeindex)

    def test_0(self):
        at = Atom()
        
        print(at)
        print(at.atomic_number)
        
        top = Topology("./data/HP36.top")
        a0 = top.atomlist[0]
        a100 = top.atomlist[100]
        a0cp = a0.copy()
        print(a0)
        print(a0cp)
        print(a100)
        assert a0.resnum == 0
        print(a100.resnum)
        print(dir(a100))
        print(a100.mass)
        print(a100.name)
        print(len(a100.name))
        print(a100.is_bonded_to(120))
        print(a100.n_bonds)

    def test_0(self):
        top = Topology("./data/HP36.top")
        for atom in top:
            print (atom)

if __name__ == "__main__":
    unittest.main()
