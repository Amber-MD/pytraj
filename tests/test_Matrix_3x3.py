import numpy as np
import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.Matrix_3x3 import Matrix_3x3

class TestMatrix_3x3(unittest.TestCase):
    def test_0(self):
        mat = Matrix_3x3(list(range(9)))
        print(mat)
        print(mat.tolist())
        print(mat[:])

        assert mat[:].shape == (3, 3)
        assert mat.buffer1d.shape == (9,)
        assert mat.tolist() == [float(x) for x in range(9)]

        mat[0, 0] = 100.
        assert mat.tolist()[0] == 100.

        npmat = np.asmatrix(mat[:])
        mat[0, 0] = 20.
        assert npmat[0, 0] == 20.
        mat *= mat
        print(mat)
        print(mat.tolist())
        mat.pprint()
        v1 = mat.row1()

        v1np = np.asarray(v1[:]).reshape(3, 1)
        print(v1np)
        print((npmat * v1np))
        print((mat * v1).tolist())
        print(mat[0, :])
        
if __name__ == "__main__":
    unittest.main()
