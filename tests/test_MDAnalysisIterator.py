from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir, duplicate_traj
import pytraj.common_actions as pyca
from pytraj.compat import zip

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        from MDAnalysis import Universe
        from pytraj.trajs.TrajectoryMDAnalysisIterator import (
                TrajectoryMDAnalysisIterator as MDIterator)

        u = Universe(traj.top.filename, traj.filename, format='mdcrd', topology_format='prmtop')
        traj_converted = mdio.load_MDAnalysis(u, top=traj.top)
        aa_eq(traj.xyz, traj_converted.xyz)
        aa_eq(traj[0].xyz, traj_converted[0].xyz)
        print (u.atoms)
        u_traj = u.trajectory

        titer = MDIterator(u, top=traj.top)
        print (titer)
        assert titer.n_atoms == u_traj.numatoms
        assert titer.n_frames == u_traj.numframes

        # make sure titer.top is Topology object
        assert isinstance(titer.top, Topology)
        print (titer[0][0])

        # make sure we get correct frame with given index
        aa_eq(traj_converted[0].xyz, titer[0].xyz)

        # test slicing
        aa_eq(titer[0:10:2].xyz, traj[0:10:2].xyz)

        # make sure we can reprodue pytraj' xyz coords
        for f0, f1 in zip(titer, traj):
            #print (f0[0])
            assert isinstance(f0, Frame) == True
            aa_eq(f0.xyz, f1.xyz)

        # let's do some analysis
        d_mda = pyca.search_hbonds(titer, dtype='ndarray')
        d_traj = pyca.search_hbonds(traj, dtype='ndarray')
        aa_eq(d_mda, d_traj)

if __name__ == "__main__":
    unittest.main()
