from __future__ import print_function
import unittest
from pytraj import io as mdio
from pytraj.compat import izip
from pytraj.core import Atom


class Test(unittest.TestCase):

    def test_0(self):
        top = mdio.load_topology("./data/Tc5b.top")

        # number
        assert isinstance(top[0], Atom)
        assert isinstance(top[:2], list)
        assert isinstance(top[:1], Atom)
        assert top[0].name == top['@1'].name

        # mask, AtomMask, python array, list
        atm = top("@CA")
        indices = atm.indices
        for a1, a2, a3, a4 in izip(top['@CA'], top[atm], top[indices], top[list(indices)]):
            assert a1.name == a2.name == a3.name == a4.name == 'CA  '

        # check len
        assert top[:].__len__() == top.n_atoms
        assert top[:10].__len__() == 10

if __name__ == "__main__":
    unittest.main()
