import unittest
from pytraj.decorators import no_test
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):
    def test_0(self):
        atom = Atom()
        print(atom.get_all_atomic_elements())
        alist = Atom.get_all_atomic_elements()
        print(alist)
        assert atom.get_bond_length('NITROGEN', 'IODINE') == 1.6
        assert atom.get_bond_length('nitrogen', 'iodine') == 1.6

    def test_build(self):
        print('test_build')
        atom = Atom("CA", "CA")
        print(atom)
        print(atom.resnum)
        print(atom.mass)
        print(atom.gb_radius)
        print(atom.element)
        print(atom.atomic_number)
        print(dir(atom))
        assert isinstance(atom, Atom) == True

        print("test copy")
        a = Atom(atom)
        print(a)

    def test_bonds(self):
        top = mdio.load_topology("./data/Tc5b.top")
        atom = top[20]
        print(atom)
        bonded_indices = atom.bonded_indices()
        for i in bonded_indices:
            assert atom.is_bonded_to(i) == True


if __name__ == "__main__":
    unittest.main()
