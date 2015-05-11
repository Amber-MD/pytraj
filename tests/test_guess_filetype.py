from __future__ import print_function
import unittest
from pytraj.testing import is_linux

class Test(unittest.TestCase):
    def test_0(self):
        if is_linux():
            from pytraj.misc import file_type_info
            fname0 = "data/Tc5b.top"
            print (file_type_info(fname0).decode())
        else:
            print ("only do this test for Linux")

if __name__ == "__main__":
    unittest.main()
