import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.math import Matrix_3x3, Vec3


class Test(unittest.TestCase):

    def test_0(self):
        mat = Matrix_3x3(list(range(10, 19)))
        vec = Vec3(list(range(3)))


if __name__ == "__main__":
    unittest.main()
