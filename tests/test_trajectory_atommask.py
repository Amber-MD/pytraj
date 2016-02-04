import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):

    def test_0(self):
        mask = "@CA"
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        top = traj.top
        atm = traj.top(mask)
        n_selected_atoms = atm.n_atoms
        newtraj = traj[atm]
        newtraj2 = traj[mask + ' :frame']
        assert_almost_equal(newtraj2.xyz.flatten(), newtraj.xyz.flatten())
        assert_almost_equal(newtraj2.xyz.flatten(), newtraj.xyz.flatten())
        assert (newtraj.xyz.shape == (traj.n_frames, n_selected_atoms, 3))

        # check if there is segmentation fault
        i = 0
        for frame in traj:
            i += 1

    def test_1(self):
        # why Trajectory is here? because I am lazy to move
        mask = "@CA"
        # creat Trajectory ( [:] )
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")[:]
        top = traj.top
        atm = traj.top(mask)
        newtraj = traj[atm]
        newtraj2 = traj[mask + ' :frame']
        assert_almost_equal(newtraj2.xyz.flatten(), newtraj.xyz.flatten())
        assert_almost_equal(newtraj2.xyz.flatten(), newtraj.xyz.flatten())


if __name__ == "__main__":
    unittest.main()
