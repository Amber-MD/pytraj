import unittest
from pytraj import *


class Test(unittest.TestCase):

    def test_0(self):
        def test_class(self):
            class FA(Trajectory):
                pass

            FA("./data/Tc5b.x", "./data/Tc5b.top")

        self.assertRaises(TypeError, lambda: test_class())


if __name__ == "__main__":
    unittest.main()
