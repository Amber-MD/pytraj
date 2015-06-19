import unittest
from pytraj.base import *
from pytraj import TrajectoryIterator
from pytraj.data_sample.load_sample_data import load_sample_data
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal as aa_eq

class Test(unittest.TestCase):
    def test_0(self):
        traj = load_sample_data()[:]
        assert isinstance(traj, Trajectory) == True
        assert traj.top.n_atoms == 34
        assert traj.shape == (1, 34, 3)

        traj2 = mdio.load_sample_data()
        assert isinstance(traj2, TrajectoryIterator) == True
        assert traj2.top.n_atoms == 34
        assert traj2.shape == (1, 34, 3)
        aa_eq(traj.xyz, traj2.xyz)

if __name__ == "__main__":
    unittest.main()
