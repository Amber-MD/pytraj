from __future__ import print_function
import unittest
from pytraj import io as mdio


class Test(unittest.TestCase):

    def test_0(self):
        framearray = mdio.iterload("./data/tz2.truncoct.nc",
                                   "./data/tz2.truncoct.parm7")[:]
        for frame in framearray:
            assert frame.has_box() == True

        framearray.set_nobox()
        for frame in framearray:
            assert frame.has_box() == False

if __name__ == "__main__":
    unittest.main()
