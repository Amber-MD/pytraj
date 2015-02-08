import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj import calculate

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/tz2.crd", "./data/tz2.parm7")[:]
        # NOT SURE CORRECTLY YET
        # FIXME
        #d0 = calculate('pairdist', 'mask "*" delta 0.1', traj)
        #print (d0)

if __name__ == "__main__":
    unittest.main()
