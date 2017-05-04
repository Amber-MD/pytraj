from __future__ import print_function
import unittest

import pytraj as pt
from pytraj.utils import Timer

from utils import fn


class Test(unittest.TestCase):
    def test_0(self):
        with Timer() as t:
            pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))


if __name__ == "__main__":
    unittest.main()
