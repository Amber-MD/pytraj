import unittest
import pytraj as pt
from pytraj.testing import aa_eq


class TestLoadSampleData(unittest.TestCase):

    def test_load_samples(self):
        traj = pt.load_sample_data()[:]
        assert isinstance(traj, pt.Trajectory) == True
        assert traj.top.n_atoms == 34
        assert traj.shape == (1, 34, 3)

        traj2 = pt.load_sample_data()
        assert isinstance(traj2, pt.TrajectoryIterator) == True
        assert traj2.top.n_atoms == 34
        assert traj2.shape == (1, 34, 3)
        aa_eq(traj.xyz, traj2.xyz)


if __name__ == "__main__":
    unittest.main()
