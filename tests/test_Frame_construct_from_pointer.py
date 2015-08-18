from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):

    def test_0(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        print(traj)

        for xyz0 in traj.xyz:
            frame = pt.Frame(traj.n_atoms, xyz0, _as_ptr=True)
            aa_eq(frame.xyz, xyz0)
            print(frame.xyz[0, 0])

    def test_1(self):
        print('api.Trajectory iterating')
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        api_t0 = pt.api.Trajectory(traj)
        print(api_t0)

        # __iter__
        for f in api_t0:
            pass

        f.xyz[0, 0] = 10.
        assert f.xyz[0, 0] == 10.
        assert api_t0.xyz[-1, 0, 0] == 10.

        # __getitem__
        # make a new copy
        api_t0 = pt.api.Trajectory(traj)
        f0 = api_t0[0]
        f0.xyz[0, 0] = 200.
        assert api_t0.xyz[0, 0, 0] == 200.

if __name__ == "__main__":
    unittest.main()
