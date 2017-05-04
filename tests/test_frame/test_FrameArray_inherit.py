import unittest
from pytraj import *


class Test(unittest.TestCase):
    def test_0(self):
        def test_class(self):
            class FA(Trajectory):
                pass

            FA(fn('Tc5b.x'), fn('Tc5b.top'))

        self.assertRaises(TypeError, lambda: test_class())


if __name__ == "__main__":
    unittest.main()
