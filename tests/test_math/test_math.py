import unittest

from pytraj.math import Matrix_3x3, Vec3


class Test(unittest.TestCase):
    def test_0(self):
        Matrix_3x3(list(range(10, 19)))
        Vec3(list(range(3)))


if __name__ == "__main__":
    unittest.main()
