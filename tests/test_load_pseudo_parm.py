from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.misc import get_atts
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import no_test, test_if_having

class Test(unittest.TestCase):
    @test_if_having('mdtraj')
    def test_0(self):
        print ("load mdtraj parm")
        import mdtraj as md
        from pytraj.externals import _load_pseudo_parm
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        m_top = md.load_prmtop("./data/Tc5b.top")
        print (m_top)
        top = mdio.load_pseudo_parm(m_top) 
        print (top)

    @test_if_having('chemistry')
    def test_1(self):
        print ("load ParmEd")
        p_top = mdio.load_ParmEd("./data/Tc5b.top")
        top = mdio.load_pseudo_parm(p_top)
        print (top)
        top.summary()

        # make sure pseudo_top is usable
        from pytraj.common_actions import calc_distance
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        fa = traj[:]
        fa.top = top.copy()

        # usingn cpptraj Topology
        d0 = calc_distance(traj, ":2@CA :10@CA")

        # using pseudo_top
        d1 = calc_distance(fa, ":2@CA :10@CA")

        assert_almost_equal(d0[:], d1[:])

if __name__ == "__main__":
    unittest.main()
