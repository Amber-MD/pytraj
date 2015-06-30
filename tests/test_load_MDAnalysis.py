from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal, eq
from pytraj.utils.check_and_assert import assert_almost_equal as aa_eq
from pytraj.decorators import no_test, test_if_having
import numpy as np


class Test(unittest.TestCase):
    #@no_test

    @test_if_having("MDAnalysis")
    def test_0(self):
        from MDAnalysis import Universe

        parm_name = "./data/ala3.psf"
        mdx_name = "./data/ala3.dcd"

        # load pytraj object
        traj = mdio.iterload(mdx_name, parm_name)
        print(traj)

        # load MDAnalysis object
        u = Universe(parm_name, mdx_name)

        # try converting MDAnalysis object to pytraj object
        md_traj = mdio.load_MDAnalysis(u)
        print(md_traj)

        # assert
        assert traj.n_frames == md_traj.n_frames
        assert traj.top.n_atoms == md_traj.top.n_atoms

        # Atom mask selection
        p_indices = traj.top("@CA").indices
        m_indices = md_traj.top("@CA").indices

        eq(p_indices, m_indices)

    #@no_test
    @test_if_having("MDAnalysis")
    @test_if_having("numpy")
    def test_1(self):
        # Aim: not getting segmentation fault
        from pytraj import io
        from MDAnalysis import Universe
        from MDAnalysisTests.datafiles import PSF, DCD
        u = Universe(PSF, DCD)
        ftraj = io.load_MDAnalysis(u)

        top = ftraj.top
        print(ftraj, top)
        print(top.box, top.box.type)
        aa_eq(sorted(top._bonds_ndarray.flatten()),
              sorted(np.array(u.universe.bonds.to_indices()).flatten()))
        print(top.residue_names)
        print(top.atom_names)
        print('n_mols = %s' % top.n_mols)
        aa_eq(sorted(top._angles_ndarray.flatten()),
              sorted(np.array(u.universe.angles.to_indices()).flatten()))
        aa_eq(sorted(top._dihedrals_ndarray.flatten()),
              sorted(np.array(u.universe.torsions.to_indices()).flatten()))

        # just want to make sure not getting segmentation fault
        ftraj.calc_COG().tolist()[:10]
        ftraj.calc_radgyr().tolist()[:10]

        # FIXME : segmentation fault, need to check cpptraj code
        d = ftraj.calc_dihedral(':2 :3 :7 :9')
        print(d.to_ndarray())

        # FIXME : segmentation fault, need to check cpptraj code
        ftraj.search_hbonds()

if __name__ == "__main__":
    unittest.main()
