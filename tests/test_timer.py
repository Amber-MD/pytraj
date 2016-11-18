from __future__ import print_function
import unittest

from pytraj import io as mdio
from pytraj.utils import Timer


class Test(unittest.TestCase):

    def test_0(self):
        with Timer() as t:
            mdio.iterload(fn('Tc5b.x'), fn('Tc5b.top'))


if __name__ == "__main__":
    unittest.main()
