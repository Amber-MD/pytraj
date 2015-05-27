from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj import calculate
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import no_test, test_if_having

class Test(unittest.TestCase):
    @no_test
    def test_0(self):
        print (calculate)
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist = calculate(traj, "distance", ":2@CA :10@CA")
        d0 = dslist[0]
        assert isinstance(d0.tolist(), list)
        assert_almost_equal(d0.data, d0.tolist())

    @test_if_having("numpy")
    #@no_test
    def test_1(self):
        print (calculate)
        import numpy as np
        print ("test dataset double")
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist = calculate("distance", traj, ":2@CA :10@CA")
        d0 = dslist[0]
        assert isinstance(d0.to_ndarray(), np.ndarray)
        assert_almost_equal(d0.data, d0.to_ndarray())

        # test memview
        arr0 = d0.to_ndarray()
        arr0[0] = 1000.
        assert arr0[0] == 1000.
        assert d0[0] == 1000.
        assert d0.data[0] == 1000.
        
    @test_if_having("numpy")
    def test_2(self):
        print (calculate)
        import numpy as np
        print ("test dataset int")
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        d0 = calculate("dssp", traj, ":2-4")['TYR:3']
        assert isinstance(d0.to_ndarray(), np.ndarray)
        assert_almost_equal(d0.data, d0.to_ndarray())

        ## test memview
        arr0 = d0.to_ndarray()
        arr0[0] = 1000
        assert arr0[0] == 1000
        assert d0[0] == 1000
        assert d0.data[0] == 1000

if __name__ == "__main__":
    unittest.main()
