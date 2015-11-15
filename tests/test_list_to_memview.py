from __future__ import print_function
import unittest
from pytraj.cyutils import _int_array1d_like_to_memview
from pytraj.cyutils import _int_array2d_like_to_memview


class Test(unittest.TestCase):

    def test_0(self):
        import numpy as np
        from array import array

        my1dlist = [100, 200]
        my1dnp = np.array(my1dlist)

        # 1D-list
        arr0 = _int_array1d_like_to_memview(my1dlist)
        atolist = list(arr0)
        assert my1dlist == atolist

        # 1D-ndarray
        arr0 = _int_array1d_like_to_memview(my1dnp)
        atolist = list(arr0)
        assert my1dnp.tolist() == atolist

        # 1D list of python array
        arr0 = _int_array1d_like_to_memview(array('i', my1dlist))
        atolist = list(arr0)
        assert my1dlist == atolist

    def test_1(self):
        import numpy as np
        from array import array

        my2dlist = [[100, 200], [1, 2]]
        my2dnp = np.array(my2dlist)

        # 2D-list
        arr0 = _int_array2d_like_to_memview(my2dlist)
        atolist = [list(x) for x in arr0]
        assert my2dlist == atolist

        # 2D-ndarray
        arr0 = _int_array2d_like_to_memview(my2dnp)
        atolist = [list(x) for x in arr0]
        assert my2dnp.tolist() == atolist

        # 2D list of python array
        arr0 = _int_array2d_like_to_memview([array('i', x) for x in my2dlist])
        atolist = [list(x) for x in arr0]
        assert my2dlist == atolist

    def test_2(self):
        # range
        arr = _int_array1d_like_to_memview(range(100))
        assert (list(arr)) == list(range(100))


if __name__ == "__main__":
    unittest.main()
