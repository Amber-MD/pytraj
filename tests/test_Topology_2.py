import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_get_iter(self):
        print("test_get_iter")
        top = Topology("./data/DOPC.parm7")
        print(top.get_resname_set())
        print(top.get_atomname_set())
        s = [atom.name for atom in top[":PC@H*"]]
        atom0 = top[":PC@H*"][0]
        print(dir(atom0))
        print(atom0.resnum)

if __name__ == "__main__":
    unittest.main()
