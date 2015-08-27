import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import test_if_having, no_test


class Test(unittest.TestCase):
    @test_if_having("numpy")
    def test_0(self):
        #print("test TrajectoryIterator")
        mask = "@CA"
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        #print(traj)
        top = traj.top
        atm = traj.top(mask)
        n_selected_atoms = atm.n_atoms
        newtraj = traj[atm]
        newtraj2 = traj[mask + ' :frame']
        #print(newtraj.xyz.shape)
        assert_almost_equal(newtraj2.xyz.flatten(), newtraj.xyz.flatten())
        assert_almost_equal(newtraj2.xyz.flatten(), newtraj.xyz.flatten())
        assert (newtraj.xyz.shape == (traj.n_frames, n_selected_atoms, 3))

        # check if there is segmentation fault
        i = 0
        for frame in traj:
            #print(frame)
            i += 1
        #print(i)

        #print(traj[:])
        #print(traj[:3])

        #@no_test
    @test_if_having("numpy")
    def test_1(self):
        # why Trajectory is here? because I am lazy to move
        #print("test Trajectory")
        mask = "@CA"
        # creat Trajectory ( [:] )
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")[:]
        #print(traj)
        top = traj.top
        atm = traj.top(mask)
        newtraj = traj[atm]
        newtraj2 = traj[mask + ' :frame']
        assert_almost_equal(newtraj2.xyz.flatten(), newtraj.xyz.flatten())
        assert_almost_equal(newtraj2.xyz.flatten(), newtraj.xyz.flatten())

        # check if there is segmentation fault
        #print(traj[0])


if __name__ == "__main__":
    unittest.main()
