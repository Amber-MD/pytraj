import unittest
from pytraj.compat import zip
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.trajs.Trajin import Trajin
from pytraj import Trajectory
import pytraj.common_actions as pyca
from pytraj.testing import aa_eq

class Test(unittest.TestCase):
    def test_0(self):
        from pytraj.TrajectoryREMDIterator import TrajectoryREMDIterator
        top = Topology("./data/Test_RemdTraj/ala2.99sb.mbondi2.parm7")

        # load regular traj
        straj = mdio.iterload("./data/Test_RemdTraj/rem.nc.000", top)

        # load all traj and extract frames having 300.0 K
        traj = mdio.load_remd("./data/Test_RemdTraj/rem.nc.000", top, "300.0")
        trajiter = mdio.iterload_remd("./data/Test_RemdTraj/rem.nc.000", top, "300.0")
        assert isinstance(traj, Trajectory)
        assert isinstance(trajiter, TrajectoryREMDIterator)
        # test slicing
        assert isinstance(trajiter[:], Trajectory)

        print (traj)
        print (trajiter)
        print (traj, traj.top, traj.n_frames)

        # make sure to get 300.0 K for all frames
        for T in traj.temperatures:
            assert_almost_equal([T], [300.0,])

        # make sure to reproduce cpptraj output
        saved_traj = mdio.iterload("data/Test_RemdTraj/temp0.crd.300.00", 
                               "./data/Test_RemdTraj/ala2.99sb.mbondi2.parm7")

        print (traj.n_frames)
        count = 0
        for f0, f1, f2 in zip(traj, trajiter, saved_traj):
            aa_eq(f0.xyz, f1.xyz)
            aa_eq(f0.xyz, f2.xyz)
            count += 1

        assert count == saved_traj.n_frames

        # test methods
        aa_eq(pyca.rmsd(trajiter), saved_traj.calc_rmsd())
        aa_eq(pyca.calc_COM(trajiter).to_ndarray(), 
              saved_traj.calc_COM().to_ndarray())
        aa_eq(pyca.calc_COG(trajiter).to_ndarray(), 
              saved_traj.calc_COG().to_ndarray())
        aa_eq(pyca.calc_dssp(trajiter, dtype='int')[0],
              saved_traj.calc_dssp(dtype='int')[0])

if __name__ == "__main__":
    unittest.main()
