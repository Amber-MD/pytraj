"""for testing bugs
"""
import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
import numpy as np

class TestBugs(unittest.TestCase):
    def test_0(self):
        trajcpp = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        farray = trajcpp[:]
        print(farray[0, 0, 0])

        # bugs: trajcpp[0, 0, 0] != farray[0, 0, 0] (must be equal)
        assert farray[0][0, 0] == trajcpp[0][0, 0]

        farray[0, 0, 0] = 100.10
        print(farray)
        print(farray[0, 0, 0])
        assert farray[0, 0, 0] == 100.10
        farray[0, 0, :] = [100.10, 101.1, 102.1]
        assert_almost_equal(farray[0, 0, :], [100.10, 101.1, 102.1])

if __name__ == "__main__":
    unittest.main()
