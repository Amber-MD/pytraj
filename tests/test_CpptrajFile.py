#import unittest
from pytraj.base import *
from pytraj.CpptrajFile import CpptrajFile


class TestCpptrajFile(unittest.TestCase):
    def test_0(self):
        # test "with" statement
        with CpptrajFile("Tc5b.crd", 'r') as cfile:
            print(cfile.is_open())
            print(dir(cfile))
            print(cfile.nextline())
            print(cfile.nextline())
            print(cfile.nextline())
            print(cfile.file_size() / 1000.)
            print(cfile.mode)


if __name__ == "__main__":
    unittest.main()
