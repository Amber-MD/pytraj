import numpy as np
import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.math import Matrix_3x3
from numpy.testing import assert_almost_equal as aa_eq_np


def npmat_fromlist(mlist):
    # convert 1D list to 2D matrix
    return np.asmatrix(np.array(mlist).reshape((3, 3)))


def eq_np_to_mat(npmat, mymat):
    assert np.any(npmat == mymat.to_ndmatrix()) == True


class TestMatrix_3x3(unittest.TestCase):

    def test_construct(self):
        from pytraj.math import Matrix_3x3 as M

        # 0
        m0 = M()
        nm0 = npmat_fromlist([0 for _ in range(9)])
        eq_np_to_mat(nm0, m0)

        # a number: 1
        m0 = M(1.)
        nm0 = npmat_fromlist([1. for _ in range(9)])
        eq_np_to_mat(nm0, m0)

        # memview
        m0 = M(range(9))
        nm0 = npmat_fromlist([float(x) for x in range(9)])
        eq_np_to_mat(nm0, m0)

        # nparray
        nm0 = npmat_fromlist([float(x) for x in range(9)])
        m0 = M(nm0)
        m1 = M(nm0.tolist())
        eq_np_to_mat(nm0, m0)
        eq_np_to_mat(nm0, m1)

    def test_math(self):
        mnp0 = np.asmatrix(np.random.rand(9).reshape(3, 3))
        mnp1 = np.asmatrix(np.random.rand(9).reshape(3, 3))
        m0 = Matrix_3x3(mnp0)
        m1 = Matrix_3x3(mnp1)
        mnp01 = mnp0 * mnp1
        m01 = m0 * m1
        aa_eq_np(m01.to_ndmatrix(), mnp01)

    def test_0(self):
        import numpy as np
        mat = Matrix_3x3(list(range(9)))

        assert mat[:].shape == (3, 3)
        assert mat.buffer1d.shape == (9, )
        asmat1 = np.asmatrix(np.arange(9).reshape((3, 3))).astype(np.float64)
        eq_np_to_mat(asmat1, mat)

        mat[0, 0] = 100.
        assert mat.tolist()[0][0] == 100.

        # test memview
        npmat = np.asmatrix(mat[:])
        mat[0, 0] = 20.
        assert npmat[0, 0] == 20.
        mat *= mat
        v1 = mat.row1
        assert v1.tolist() == list(mat[0])
        assert np.any(npmat == mat.to_ndmatrix()) == True

        mat_as_ndmatrix = mat.as_ndmatrix()
        mat_as_ndmatrix[0] = 10000.
        assert np.any(npmat == mat.to_ndmatrix()) == True
        assert np.any(npmat == mat_as_ndmatrix) == True

        v1np = np.asarray(v1[:]).reshape(3, 1)


if __name__ == "__main__":
    unittest.main()
