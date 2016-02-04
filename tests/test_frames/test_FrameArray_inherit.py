import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):

    def test_0(self):
        def test_class(self):
            class FA(Trajectory):
                pass

            fa = FA("./data/Tc5b.x", "./data/Tc5b.top")

        self.assertRaises(TypeError, lambda: test_class())


if __name__ == "__main__":
    unittest.main()
