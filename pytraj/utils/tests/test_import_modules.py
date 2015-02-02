import unittest
from pycpptraj.base import *
from pycpptraj import io as mdio
from pycpptraj.utils.check_and_assert import _import

class Test(unittest.TestCase):
    def test_0(self):
        has_np, np = _import('numpy')
        print has_np, np
        has_h5py, h5py = _import('h5py')
        

if __name__ == "__main__":
    unittest.main()
