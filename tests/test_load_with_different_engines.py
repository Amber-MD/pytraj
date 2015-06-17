from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir, duplicate_traj
import pytraj.common_actions as pyca

class Test(unittest.TestCase):
    @test_if_having("MDAnalysis")
    def test_0(self):
        from MDAnalysisTests.datafiles import PSF, DCD
        from pytraj import TrajectoryIterator
        from pytraj.trajs.TrajectoryMDAnalysisIterator import TrajectoryMDAnalysisIterator
        traj = mdio.iterload(DCD, PSF, engine='pytraj')
        traj2 = mdio.iterload(DCD, PSF, engine='mdanalysis')
        assert isinstance(traj, TrajectoryIterator) == True
        assert isinstance(traj2, TrajectoryMDAnalysisIterator) == True
        aa_eq(traj.xyz, traj2.xyz)

        print (pyca.calc_COM(traj2).to_ndarray())

    @test_if_having("MDAnalysis")
    @test_if_having("mdtraj")
    def test_0(self):
        from MDAnalysisTests.datafiles import PSF, DCD
        from pytraj import Trajectory
        from pytraj import TrajectoryIterator
        from pytraj.trajs.TrajectoryMDAnalysisIterator import TrajectoryMDAnalysisIterator

        filename = DCD
        topname = PSF
        traj = mdio.iterload(filename, topname, engine='pytraj')
        traj2 = mdio.load(filename, topname, engine='mdtraj')

        # test autoconvert
        mdio.load(filename, topname, engine='mdtraj', autoconvert=False)
        assert isinstance(traj, TrajectoryIterator) == True
        assert isinstance(traj2, Trajectory) == True
        aa_eq(traj.xyz, traj2.xyz)

        print (pyca.calc_COM(traj2).to_ndarray())

        # raise if use engine='mdtraj' for iterload
        self.assertRaises(ValueError, lambda : mdio.iterload(filename, topname, engine='mdtraj'))

if __name__ == "__main__":
    unittest.main()
