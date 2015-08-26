import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.testing import no_test


class Test(unittest.TestCase):
    def test_1(self):
        top = Topology("./data/HP36.top")
        a0 = top.atomlist[0]
        a100 = top.atomlist[100]
        a0cp = a0.copy()
        #print(a0)
        #print(a0cp)
        #print(a100)
        assert a0.resnum == 0
        #print(a100.resnum)
        #print(dir(a100))
        #print(a100.mass)
        #print(a100.name)
        #print(len(a100.name))
        #print(a100.is_bonded_to(120))
        #print(a100.n_bonds)

    def test_2(self):
        top = Topology("./data/HP36.top")
        for atom in top:
            #print(atom)

    def test_3(self):
        # test Atom() from `pytraj` namespace
        # Aim: test if having segmentation fault
        import pytraj as pt
        pt.Atom()
        # FIXME: segmentation fault
        #print(pt.Atom())


if __name__ == "__main__":
    unittest.main()
