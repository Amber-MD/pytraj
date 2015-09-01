import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.common_actions import *


class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")[:]
        for _ in range(2):
            traj.join(traj[:])
        #print(calc_dihedral(mask=":2@CA :3@CA :10@CA :11@CA", traj=traj)[:])
        #print(calc_distance(mask=":2@CA :3@CA", traj=traj)[:])
        #print(calc_angle(mask=":2@CA :3@CA :10@CA", traj=traj)[:])
        #print(calc_radgyr(mask="@CA", traj=traj)[:])
        #print(calc_molsurf(mask="@CA", traj=traj)[:])
        #print(type(calc_molsurf(mask="@CA", traj=traj)[:]))

    def test_1(self):
        #print("test mix traj/frame")
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")[:]
        d0 = calc_molsurf(
            mask="@CA",
            traj=(
                traj, traj[:2], traj[:], traj.iterframe(), traj[0]),
            dtype='dataset')
        assert d0[0].size == 33


if __name__ == "__main__":
    unittest.main()
