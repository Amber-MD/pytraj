import unittest
import pytraj as pt
from pytraj.base import *
from pytraj.utils.check_and_assert import assert_almost_equal, eq


class TestAtomMask(unittest.TestCase):

    def test_0(self):
        atm = AtomMask("@CA")
        assert atm.n_atoms == 0
        top = pt.load_topology("./data/Tc5b.top")
        top.set_integer_mask(atm)
        assert atm.n_atoms == 20
        top2 = top._modify_state_by_mask(atm)
        assert top2.n_atoms == 20
        for atom in top2:
            assert atom.name == 'CA'

        atm.invert_mask()
        top3 = top._modify_state_by_mask(atm)
        assert top3.n_atoms == top.n_atoms - 20

    def test_1(self):
        atm = AtomMask(10)
        assert atm.n_atoms == 1

    def test_3_indexing(self):
        top = pt.load_topology("./data/Tc5b.top")
        atm = AtomMask("@CA")
        top.set_integer_mask(atm)

    def test_4(self):
        from array import array
        atm = AtomMask()
        indices = array('i', range(100))
        atm.add_selected_indices(indices)
        assert_almost_equal(indices, atm.selected_indices())
        assert_almost_equal(indices, atm.indices)

        # test list
        atm2 = AtomMask(list(indices))
        assert_almost_equal(indices, atm2.indices)

        # test range
        r100 = range(100)
        atm3 = AtomMask(range(100))
        assert_almost_equal(indices, atm3.indices)

    def test_5_speed(self):
        top = pt.load_topology("./data/DOPC.parm7")
        from time import time

        t0 = time()
        indices = top(":WAT").indices
        gap_0 = time() - t0

        t0 = time()
        _indices_view = top(":WAT")._indices_view
        gap_1 = time() - t0

        count = 0
        for i, j in zip(indices, _indices_view):
            if not i == j:
                count += 1
        assert count == 0

    def test_6_speed(self):
        import numpy as np
        from pytraj import AtomMask
        # test constructor from list/array/python array
        top = pt.load_topology("./data/DOPC.parm7")
        indices = top.select(":WAT")

        atm1 = AtomMask(indices)
        atm2 = AtomMask(list(indices))
        atm3 = AtomMask(np.asarray(indices))
        # use max_atoms
        atm4 = AtomMask(np.asarray(indices), 1000)

        import numpy as np
        assert np.all(atm1.indices == atm2.indices)
        assert np.all(atm3.indices == atm4.indices)
        # FIXME: can not catch RuntimeError here
        # since we don't set atm3 max_atoms, we expect to get RuntimeError
        # if using invert_mask
        # TODO: assert fails


if __name__ == "__main__":
    unittest.main()
