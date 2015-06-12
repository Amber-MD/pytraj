from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.misc import get_atts
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.testing import aa_eq
from pytraj.decorators import no_test, test_if_having
from pytraj.compat import zip

class Test(unittest.TestCase):
    @test_if_having('mdtraj')
    def test_0(self):
        print ("load mdtraj parm")
        import mdtraj as md
        from pytraj.externals import _load_pseudo_parm
        traj = mdio.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        m_top = md.load_prmtop("./data/tz2.ortho.parm7")
        print (m_top)
        top = mdio.load_pseudo_parm(m_top) 
        print (top)
        aa_eq(top.mass, traj.top.mass, decimal=2)

    @test_if_having('parmed')
    def test_1(self):
        print ("load ParmEd")
        p_top = mdio._load_parmed("./data/Tc5b.top")
        top = mdio.load_pseudo_parm(p_top)
        print (top)
        top.summary()
        # make sure pseudo_top is usable
        from pytraj.common_actions import calc_distance
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        fa = traj[:]
        fa.top = top.copy()

        # usingn cpptraj Topology
        d0 = calc_distance(traj, ":2@CA :10@CA")

        # using pseudo_top
        d1 = calc_distance(fa, ":2@CA :10@CA")

        assert_almost_equal(d0[:], d1[:])
        aa_eq(top.mass, traj.top.mass)
        for a0, a1 in zip(top.atoms, traj.top.atoms):
            assert a0.mass == a1.mass
            assert abs(a0.charge - a1.charge) < 1E-3

    @test_if_having('parmed')
    def test_2(self):
        from pytraj import io
        import parmed as pmd
        print (pmd)
        pdb = pmd.download_PDB("1o15")
        top = io.load_pseudo_parm(pdb, guess_bond=False)
        ptop = io.loadpdb_rcsb("1o15").top


if __name__ == "__main__":
    unittest.main()
