from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.testing import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca

class Test(unittest.TestCase):
    @test_if_having("numpy")
    def test_0(self):
        import numpy as np
        parm_name = "./data/Tc5b.top"
        true_top = mdio.load(parm_name)

        # load ParmEd
        parm = mdio.load_ParmEd(parm_name) 

        # load pseudo_parm
        ptop = mdio.load_pseudo_parm(parm)
        print (ptop)

        # assert
        eq(sorted(ptop._bonds_ndarray.flatten()), 
           sorted(true_top._bonds_ndarray.flatten()))

        eq(sorted(ptop._angles_ndarray.flatten()), 
           sorted(true_top._angles_ndarray.flatten()))

        eq(sorted(ptop._dihedrals_ndarray.flatten()), 
           sorted(true_top._dihedrals_ndarray.flatten()))

        assert (ptop.box.type == 'nobox')

if __name__ == "__main__":
    unittest.main()
