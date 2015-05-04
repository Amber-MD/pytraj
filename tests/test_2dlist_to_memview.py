from __future__ import print_function
import unittest

class Test(unittest.TestCase):
    def test_0(self):
        from pytraj._utils import _list2d_to_memview
        mylist = [[100, 200], [1, 2]]
        arr0 = _list2d_to_memview(mylist)
        assert hasattr(arr0, 'memview') == True
        atolist = [list(x) for x in arr0]
        assert mylist == atolist

if __name__ == "__main__":
    unittest.main()
