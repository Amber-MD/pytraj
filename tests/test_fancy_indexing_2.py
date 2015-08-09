import unittest
from pytraj import io as mdio
from pytraj import Trajectory
from pytraj.utils.check_and_assert import assert_almost_equal as aa_eq


class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        FA = traj[:]

        print(traj['@CA'])
        frame0 = traj[0]
        print(hasattr(frame0, 'shape'))
        aa_eq(frame0[traj.top("@CA")].flatten(), traj['@CA'].xyz.flatten())

        # slicing with list or array
        indices = [1, 2, 3]
        fa = traj[indices]
        fa2 = FA[indices]
        fa3 = traj[range(1, 4)]
        fa4 = FA[range(1, 4)]
        self.assertIsInstance(fa, Trajectory)
        # from TrajectoryIterator
        aa_eq(fa[0].coords, traj[1].coords)
        aa_eq(fa[1].coords, traj[2].coords)
        # from Trajectory
        aa_eq(fa2[1].coords, traj[2].coords)
        aa_eq(fa2[0].coords, traj[1].coords)

        # from "range"
        aa_eq(fa3[1].coords, traj[2].coords)
        aa_eq(fa3[0].coords, traj[1].coords)
        aa_eq(fa4[1].coords, traj[2].coords)
        aa_eq(fa4[0].coords, traj[1].coords)

    def test_1(self):
        # AtomMask
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        fa = traj.to_mutable_trajectory()
        xyz = traj.xyz[:]
        atm = traj.top.select("@CA")
        indices = atm.indices

        aa_eq(fa[0, atm], fa[0][atm])
        aa_eq(traj[0, atm], fa[0][atm])
        aa_eq(traj[0, atm, 0], fa[0][atm, 0])
        aa_eq(traj[0, atm, 0], xyz[0][indices][0])
        aa_eq(traj[0, '@CA', 0], xyz[0][indices][0])


if __name__ == "__main__":
    unittest.main()
