from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir

class Test(unittest.TestCase):
    @test_if_having("mdtraj")
    def test_0(self):
        import pytraj.common_actions as ca
        import mdtraj as md
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        m_top = md.load_prmtop("./data/Tc5b.top")
        m_traj = md.load_mdcrd("./data/md1_prod.Tc5b.x", m_top)
        v0 = ca.calc_vector("@CA @N,C,O mass", traj)

if __name__ == "__main__":
    unittest.main()
