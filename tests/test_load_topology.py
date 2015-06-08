from __future__ import print_function
import unittest
from pytraj import io as mdio

class Test(unittest.TestCase):
    def test_0(self):
        top0 = mdio.load("./data/Tc5b.top")
        top1 = mdio.load_topology(top0.filename)
        assert top0.n_atoms == top1.n_atoms

if __name__ == "__main__":
    unittest.main()
