from __future__ import print_function
import unittest
from pytraj import io as mdio
from pytraj.utils import aa_eq
from pytraj.decorators import no_test


class Test(unittest.TestCase):
    #@no_test

    def test_0(self):
        from pytraj import Trajectory
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        trajiter = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        print(traj)
        t = Trajectory()
        t.top = traj.top.copy()
        t.append(traj[0], copy=False)
        print(t)
        #frame = t[0]
        frame = t[0]
        print(frame.n_atoms)
        print(frame[0])

        my_number = [1000. for _ in range(3)]
        t[0, 0] = my_number
        aa_eq(t[0, 0], my_number)
        aa_eq(traj[0, 0], my_number)
        print(traj[0, 0], t[0, 0])

    def test_1(self):
        from pytraj import Trajectory
        itraj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        t = Trajectory()
        t.top = itraj.top.copy()
        t.append(itraj[0])
        aa_eq(t.xyz, itraj.xyz)

    def test_2(self):
        from pytraj import Trajectory
        itraj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        t = Trajectory()
        t.top = itraj.top.copy()
        t.append(itraj[0].copy())
        aa_eq(t.xyz, itraj.xyz)


if __name__ == "__main__":
    unittest.main()
