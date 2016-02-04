"""for testing bugs
"""
import unittest
import numpy as np
import pytraj as pt
from pytraj.testing import aa_eq


class TestBugs(unittest.TestCase):

    def test_0(self):
        trajcpp = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        farray = trajcpp[:]

        # bugs: trajcpp[0, 0, 0] != farray[0, 0, 0] (must be equal)
        assert farray[0][0, 0] == trajcpp[0][0, 0]

        farray[0, 0, 0] = 100.10
        assert farray[0, 0, 0] == 100.10
        farray[0, 0, :] = [100.10, 101.1, 102.1]
        aa_eq(farray[0, 0, :], [100.10, 101.1, 102.1])


if __name__ == "__main__":
    unittest.main()
