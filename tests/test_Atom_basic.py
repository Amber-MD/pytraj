from __future__ import print_function
import unittest
import pytraj as pt


class Test(unittest.TestCase):

    def test_0(self):
        # check segmentation fault
        at = pt.Atom()


if __name__ == "__main__":
    unittest.main()
