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
    @test_if_having("chemistry")
    def test_0(self):
        import numpy as np
        import chemistry as chem
        parm_name = "./data/Tc5b.top"
        traj = mdio.iterload("./data/md1_prod.Tc5b.x",  parm_name)
        true_top = mdio.iterload(parm_name)

        # load ParmEd
        parm = mdio._load_chem(parm_name) 
        assert isinstance(parm, chem.Structure)
        parm.load_coordinates(traj[0].coords)

        # load pseudo_parm
        ptop = mdio.load_pseudo_parm(parm)
        fake_fa = mdio.load_ParmEd(parm, restype='traj')
        assert isinstance(fake_fa, Trajectory)
        aa_eq(fake_fa[0].coords, traj[0].coords)
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
