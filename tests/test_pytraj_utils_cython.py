import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj._utils import test
import numpy as np

class Test(unittest.TestCase):
    def test_0(self):
        # FIXME: wrong result. got 0.0 for all elements.
        # I expect got memview array of vector[double] = range(100)
        print (test)
        arr0 = (test())
        print (np.asarray(arr0))

if __name__ == "__main__":
    unittest.main()
