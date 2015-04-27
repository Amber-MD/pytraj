from __future__ import print_function
import unittest
from pytraj import io as mdio
from pytraj.testing import test_if_having
from pytraj.utils import assert_almost_equal as aa_e

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        top = traj.top

        # numbers wer taken from top.summary() (C++ stdout)
        assert list(top.bonds).__len__() == 310 
        assert list(top.angles).__len__() == 565 
        assert list(top.dihedrals).__len__() == 1351

    @test_if_having("numpy")
    @test_if_having("chemistry")
    def test_1(self):
        # try to rebuild pseudoparm
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        top = traj.top

        # load ParmEd object
        parm = mdio.load_ParmEd("./data/Tc5b.top")

        # make pseudo parm
        ptop = mdio.load_pseudo_parm(parm)
        ptop.add_bonds(top._bonds_ndarray)
        aa_e(ptop._bonds_ndarray.ravel(), top._bonds_ndarray.ravel())

        ptop.add_angles(top._angles_ndarray)
        aa_e(ptop._angles_ndarray.ravel(), top._angles_ndarray.ravel())

        ptop.add_dihedrals(top._dihedrals_ndarray)
        aa_e(ptop._dihedrals_ndarray.ravel(), top._dihedrals_ndarray.ravel())


if __name__ == "__main__":
    unittest.main()
