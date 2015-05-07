import unittest
from pytraj.base import *
from pytraj.decorators import no_test
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    #@no_test
    def test_0(self):
        traj = mdio.iterload("./data/tz2.truncoct.nc", "./data/tz2.truncoct.parm7")
        act = adict['radial']
        dslist = DataSetList()
        # why need `radial` keyword here?
        act("radial 0.5 10.0 :5@CD :WAT@O", traj(0, 9), traj.top, dslist=dslist)
        act.print_output()
        print (dslist.size)

        print (dslist[0].size)
        print (dslist[0][:])

    def test_1(self):
        from pytraj.common_actions import calc_radial
        traj = mdio.iterload("./data/tz2.truncoct.nc", "./data/tz2.truncoct.parm7")
        d0 = calc_radial(traj(0, 9), "0.5 10.0 :5@CD :WAT@O", top=traj.top)
        print (d0)
        print (d0.size)
        print (d0[0][:])

        # assert
        try:
            import os
            import numpy as np
            amberhome =  os.environ['AMBERHOME']
            saved_data_dir = amberhome + "/AmberTools/test/cpptraj/Test_Radial/cRadial.agr.save"
            print (saved_data_dir)
            data = np.loadtxt(saved_data_dir, skiprows=8).transpose()[1]
            print (data)
            assert_almost_equal(d0[0][:], data)
        except (EnvironmentError, KeyError):
            print ("can not find cpptraj test. Skip")

if __name__ == "__main__":
    unittest.main()
