from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal, eq
from pytraj.decorators import no_test, test_if_having

class Test(unittest.TestCase):
    @test_if_having("MDAnalysis")
    def test_0(self):
        from pytraj._load_MDAnalysis import load_MDAnalysis
        from MDAnalysis import Universe

        parm_name = "./data/ala3.psf"
        mdx_name = "./data/ala3.dcd"

        # load pytraj object
        traj = mdio.load(mdx_name, parm_name)
        print (traj)

        # load MDAnalysis object
        u = Universe(parm_name, mdx_name)

        # try converting MDAnalysis object to pytraj object
        md_traj = load_MDAnalysis(u)
        print (md_traj)

        # assert
        assert traj.n_frames == md_traj.n_frames
        assert traj.top.n_atoms == md_traj.top.n_atoms

        # Atom mask selection
        p_indices = traj.top("@CA").indices
        m_indices = md_traj.top("@CA").indices 

        eq(p_indices, m_indices)

if __name__ == "__main__":
    unittest.main()
