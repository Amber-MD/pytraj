import unittest
import pytraj as pt
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):
    def test_get_iter(self):
        #print("test_get_iter")
        top = pt.load_topology("./data/DOPC.parm7")
        s = [atom.name for atom in top[":PC@H*"]]
        atom0 = top[":PC@H*"][0]

        old_natoms = top.n_atoms
        print('test')
        self.assertRaises(ValueError, lambda: top.join(top))
        top.join(top.copy())
        assert top.n_atoms == 2 * old_natoms


if __name__ == "__main__":
    unittest.main()
