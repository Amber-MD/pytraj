import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        from pytraj import datasets
        ddict = datasets.__dict__
        keys = [key for key in ddict if 'DataSet' in key]
        print (keys)

        # remove base classes
        useless_keys = ['DataSet', 'DataSet_1D', 'DataSet_2D', 'DataSet_Coords', 'DataSet_Modes']
        for _key in useless_keys:
            keys.remove(_key)

        print (keys)

        for key in keys:
            d0 = ddict[key]()
            print (d0.name, d0.dtype)

if __name__ == "__main__":
    unittest.main()
