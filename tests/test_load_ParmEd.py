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
    @test_if_having("parmed")
    def test_0(self):
        import numpy as np
        import parmed as pmd
        parm_name = "./data/Tc5b.top"
        traj = mdio.iterload("./data/md1_prod.Tc5b.x",  parm_name)
        print(traj[0].coords[:10])
        true_top = mdio.load(parm_name)

        # load ParmEd
        parm = mdio._load_parmed(parm_name)
        assert isinstance(parm, pmd.Structure)

        # load pseudo_parm
        ptop = mdio.load_pseudo_parm(parm)

        # assert
        eq(sorted(ptop._bonds_ndarray.flatten()),
           sorted(true_top._bonds_ndarray.flatten()))

        eq(sorted(ptop._angles_ndarray.flatten()),
           sorted(true_top._angles_ndarray.flatten()))

        eq(sorted(ptop._dihedrals_ndarray.flatten()),
           sorted(true_top._dihedrals_ndarray.flatten()))

        assert (ptop.box.type == 'nobox')

    @test_if_having("numpy")
    @test_if_having("parmed")
    def test_1(self):
        import numpy as np
        import parmed as pmd
        from pytraj.externals._load_ParmEd import to_ParmEd
        parm_name = "./data/Tc5b.top"
        traj = mdio.iterload("./data/md1_prod.Tc5b.x",  parm_name)
        parm = to_ParmEd(traj.top)
        top2 = mdio.load_pseudo_parm(parm)

        true_top = traj.top
        ptop = top2
        eq(sorted(ptop._bonds_ndarray.flatten()),
           sorted(true_top._bonds_ndarray.flatten()))

        eq(sorted(ptop._angles_ndarray.flatten()),
           sorted(true_top._angles_ndarray.flatten()))

        eq(sorted(ptop._dihedrals_ndarray.flatten()),
           sorted(true_top._dihedrals_ndarray.flatten()))

    @no_test
    @test_if_having("numpy")
    @test_if_having("parmed")
    def test_2(self):
        # turn off test to check loading code
        import pytraj.io as io
        # try loading PSF and doing analysis
        import numpy as np
        import parmed as pmd
        parm_name = "./data/ala3.psf"
        traj = mdio.iterload("./data/ala3.dcd",  parm_name)
        parm = pmd.load_file(parm_name)
        #p_top = io.load_pseudo_parm(parm)
        p_top = io.load_full_ParmEd(parm)
        print('test2: parm', parm.__repr__())
        print('test2: p_top', p_top)
        traj_new_ptop = mdio.iterload(traj.filename, top=p_top)
        # use `search_hbonds` since I got segfault with MDAnalysis
        ds = traj.search_hbonds(dtype='ndarray')
        ds_newtop = traj_new_ptop.search_hbonds(dtype='ndarray')
        aa_eq(ds, ds_newtop)

if __name__ == "__main__":
    unittest.main()
