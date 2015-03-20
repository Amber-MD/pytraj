from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        act = adict['multidihedral']
        dslist = DataSetList()
        dflist = DataFileList()
        act("multidihedral phi psi resrange 6-9 out ./output/test_muldih.dat", traj[:2], 
            dslist=dslist, dflist=dflist)
        act.print_output()

        print (dslist.size)
        for d0 in dslist:
            print ("d0.size, d0.dtype, d0.name: ", d0.size, d0.dtype, d0.name)
            print ("d0.scalar_type, d0.scalar_mode: ", d0.scalar_type, d0.scalar_mode)
            print ("d0.aspect, d0.legend: ", d0.aspect, d0.legend)
            print ("d0.idx, d0.name, d0.ndim: ", d0.idx, d0.name, d0.ndim)
            d0.info()

        print (dir(d0))
        print (d0.is_torsion_array())

        print (dslist.get_dtypes())
        print (dslist.get_aspects())
        print (dslist.get_scalar_modes())
        print (dslist.get_scalar_types())
        print (dslist.get_legends())

if __name__ == "__main__":
    unittest.main()
