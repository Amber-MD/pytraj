import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.NameType import NameType

class TestNameType(unittest.TestCase):
    def test_0(self):
        nt = NameType("CA")
        print(dir(nt))
        print(nt[2])

        nt2 = nt.copy()
        assert (nt2 == nt) == True
        assert nt.match(nt2) == True
        assert nt2 == 'CA'
        print(nt.truncated())

        assert nt == 'CA'

if __name__ == "__main__":
    unittest.main()
