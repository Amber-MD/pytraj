from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.externals.netcdf import netcdf_file

class Test(unittest.TestCase):
    def test_0(self):
        fh = netcdf_file("./data/DPDP.nc")
        print (fh)
        print (dir(fh))
        print (fh.variables)

if __name__ == "__main__":
    unittest.main()
