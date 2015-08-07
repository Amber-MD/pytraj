from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import has_amberhome, amberhome, cpptraj_test_dir


class Test(unittest.TestCase):
    def test_0(self):
        print("has_amberhome = ", has_amberhome)
        print("amberhome = ", amberhome)
        print("cpptraj_test_dir= ", cpptraj_test_dir)


if __name__ == "__main__":
    unittest.main()
